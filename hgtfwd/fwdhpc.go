package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/covs"
	"github.com/mingzhi/hgt/fwd"
	"log"
	"math/rand"
	"os"
	"runtime"
)

var (
	size     int     // popoulation size
	reps     int     // replications
	lens     int     // genome lengths
	maxl     int     // max l
	frag     int     // fragment length
	gens     int     // number of generations
	mutation float64 // mutation rate
	transfer float64 // transfer rate
	prefix   string  // prefix
)

type Moment struct {
	Mean *desc.Mean
	Sd   *desc.StandardDeviation
}

func (m *Moment) Increment(d float64) {
	m.Mean.Increment(d)
	m.Sd.Increment(d)
}

type Result struct {
	ks, vd                             float64
	scovs, rcovs, xyPL, xsysPL, smXYPL []float64
}

func init() {
	flag.IntVar(&size, "size", 1000, "population size")
	flag.IntVar(&lens, "genome", 1000, "genome length")
	flag.IntVar(&frag, "frag", 100, "fragment length (ratio)")
	flag.IntVar(&reps, "reps", 1000, "repeats")
	flag.IntVar(&maxl, "maxl", 100, "maxl")
	flag.IntVar(&gens, "gens", 10000, "number of generations")
	flag.Float64Var(&transfer, "transfer", 1e-4, "transfer rate")
	flag.Float64Var(&mutation, "mutation", 1e-4, "mutation rate")
	flag.StringVar(&prefix, "prefix", "test", "prefix")

	flag.Parse()
	if maxl < 2*frag {
		maxl = 2 * frag
	}
}

func main() {
	ncpu := runtime.NumCPU()
	runtime.GOMAXPROCS(ncpu)

	ch := make(chan Result, ncpu)

	for i := 0; i < ncpu; i++ {
		b := i * reps / ncpu
		e := (i + 1) * reps / ncpu
		go simulateSome(b, e, ch)
	}

	dfile, err := os.Create(fmt.Sprintf("%s_d.cov", prefix))
	if err != nil {
		log.Panic(err)
	}
	defer dfile.Close()

	dfile.WriteString(fmt.Sprintf("#size: %d\n", size))
	dfile.WriteString(fmt.Sprintf("#length: %d\n", lens))
	dfile.WriteString(fmt.Sprintf("#fragment: %d\n", frag))
	dfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
	dfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
	dfile.WriteString(fmt.Sprintf("#generations: %d\n", gens))
	dfile.WriteString("#ks, vd\n")

	moments := make([][]Moment, 5)
	for i := 0; i < len(moments); i++ {
		for j := 0; j < maxl; j++ {
			m := Moment{
				Mean: desc.NewMean(),
				Sd:   desc.NewStandardDeviationWithBiasCorrection(),
			}
			moments[i] = append(moments[i], m)
		}
	}

	for i := 0; i < reps; i++ {
		result := <-ch

		dfile.WriteString(fmt.Sprintf("%g, %g\n", result.ks, result.vd))
		for j := 0; j < maxl; j++ {
			moments[0][j].Increment(result.scovs[j])
			moments[1][j].Increment(result.rcovs[j])
			moments[2][j].Increment(result.xyPL[j])
			moments[3][j].Increment(result.xsysPL[j])
			moments[4][j].Increment(result.smXYPL[j])
		}

		if (i+1)%(reps/100) == 0 {
			cfile, err := os.Create(fmt.Sprintf("%s_covs.csv", prefix))
			if err != nil {
				log.Panic(err)
			}

			cfile.WriteString(fmt.Sprintf("#size: %d\n", size))
			cfile.WriteString(fmt.Sprintf("#length: %d\n", lens))
			cfile.WriteString(fmt.Sprintf("#fragment: %d\n", frag))
			cfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
			cfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
			cfile.WriteString(fmt.Sprintf("#generations: %d\n", gens))
			cfile.WriteString(fmt.Sprintf("#replicates: %d\n", i+1))
			cfile.WriteString("#dist, scov, rcov, xy, xsys, smxy_sd, scov_sd, rcov_sd, xy_sd, xsys_sd, smxy_sd\n")

			for j := 0; j < maxl; j++ {
				cfile.WriteString(fmt.Sprintf("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
					j,
					moments[0][j].Mean.GetResult(),
					moments[1][j].Mean.GetResult(),
					moments[2][j].Mean.GetResult(),
					moments[3][j].Mean.GetResult(),
					moments[4][j].Mean.GetResult(),
					moments[0][j].Sd.GetResult(),
					moments[1][j].Sd.GetResult(),
					moments[2][j].Sd.GetResult(),
					moments[3][j].Sd.GetResult(),
					moments[4][j].Sd.GetResult(),
				))
			}
			cfile.Close()
			log.Printf("Finish %%%d\n", (i+1)/(reps/100))
		}
	}
}

func simulateSome(b, e int, ch chan Result) {
	for i := b; i < e; i++ {
		sp := fwd.NewSeqPop(size, lens, mutation, transfer, frag)
		sp.Seed(i)
		for j := 0; j < gens; j++ {
			sp.Evolve()
		}

		seqs := sp.GetGenomes()

		diffmatrix := [][]int{}
		sample := 1000
		for j := 0; j < sample; j++ {
			a := rand.Intn(size)
			b := rand.Intn(size)
			for a == b {
				b = rand.Intn(size)
			}

			diff := []int{}

			for k := 0; k < lens; k++ {
				if seqs[a][k] != seqs[b][k] {
					diff = append(diff, k)
				}
			}

			diffmatrix = append(diffmatrix, diff)
		}

		cmatrix := covs.NewCMatrix(sample, lens, diffmatrix)

		ks, vd := cmatrix.D()
		scovs, rcovs, xyPL, xsysPL, smXYPL := cmatrix.CovCircle(maxl)
		result := Result{
			ks:     ks,
			vd:     vd,
			scovs:  scovs,
			rcovs:  rcovs,
			xyPL:   xyPL,
			xsysPL: xsysPL,
			smXYPL: smXYPL,
		}

		ch <- result
	}
}
