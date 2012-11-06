package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/coals"
	"github.com/mingzhi/hgt/covs"
	"os"
	"runtime"
	"time"
)

var (
	size     int     // population size
	sample   int     // sample size
	length   int     // genome length
	repeats  int     // number of repeats
	maxl     int     // max l
	mutation float64 // mutation rate
	transfer float64 // transfer rate
	fragment int     // transfer fragment
	prefix   string  // prefix
)

type Moments struct {
	Mean *desc.Mean
	Sd   *desc.StandardDeviation
}

func (m *Moments) Increment(d float64) {
	m.Mean.Increment(d)
	m.Sd.Increment(d)
}

type Results struct {
	ks, vd                             float64
	scovs, rcovs, xyPL, xsysPL, smXYPL []float64
}

func init() {
	flag.IntVar(&size, "size", 1000000, "population size")
	flag.IntVar(&sample, "sample", 2, "sample size")
	flag.IntVar(&length, "genome", 10000, "genome length")
	flag.IntVar(&fragment, "frag", 100, "fragment length (ratio)")
	flag.IntVar(&repeats, "rep", 1000, "repeats")
	flag.IntVar(&maxl, "maxl", 100, "maxl")
	flag.Float64Var(&transfer, "transfer", 1e-6, "transfer rate")
	flag.Float64Var(&mutation, "mutation", 1e-8, "mutation rate")
	flag.StringVar(&prefix, "prefix", "", "prefix")

	flag.Parse()
	if maxl < 2*fragment {
		maxl = 2 * fragment
	}
}

func main() {
	ncpu := runtime.NumCPU()
	runtime.GOMAXPROCS(ncpu)
	ch := make(chan Results, ncpu)

	t0 := time.Now()

	for i := 0; i < ncpu; i++ {
		begin := i * repeats / ncpu
		end := (i + 1) * repeats / ncpu
		go simusome(begin, end, ch)
	}

	analysis(ch)

	t1 := time.Now()
	fmt.Println(t1.Sub(t0))
}

func simusome(begin, end int, ch chan Results) {
	for i := begin; i < end; i++ {
		w := coals.NewWFPopulation(size, sample, length, mutation, transfer, fragment)
		w.Seed(i)
		w.Backtrace()
		seqs := w.Fortrace()
		diffmatrix := [][]int{}
		sample := w.SampleSize
		for i := 0; i < sample; i++ {
			for j := i + 1; j < sample; j++ {
				diff := []int{}
				for k := 0; k < length; k++ {
					if seqs[i][k] != seqs[j][k] {
						diff = append(diff, k)
					}
				}
				diffmatrix = append(diffmatrix, diff)
			}
		}

		cmatrix := covs.NewCMatrix(len(diffmatrix), length, diffmatrix)
		ks, vd := cmatrix.D()
		scovs, rcovs, xyPL, xsysPL, smXYPL := cmatrix.CovCircle(maxl)
		results := Results{
			ks:     ks,
			vd:     vd,
			scovs:  scovs,
			rcovs:  rcovs,
			xyPL:   xyPL,
			xsysPL: xsysPL,
			smXYPL: smXYPL,
		}
		ch <- results
	}
}

func analysis(ch chan Results) {
	dfile, err := os.Create(fmt.Sprintf("%s_d.csv", prefix))
	if err != nil {
		panic(err)
	}
	defer dfile.Close()

	dfile.WriteString(fmt.Sprintf("#size: %d\n", size))
	dfile.WriteString(fmt.Sprintf("#length: %d\n", length))
	dfile.WriteString(fmt.Sprintf("#fragment: %d\n", fragment))
	dfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
	dfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
	dfile.WriteString(fmt.Sprintf("#sample: %d\n", sample))
	dfile.WriteString("#ks, vd\n")

	momentArr := make([][]Moments, 5)
	for i := 0; i < len(momentArr); i++ {
		for j := 0; j < maxl; j++ {
			moments := Moments{
				Mean: desc.NewMean(),
				Sd:   desc.NewStandardDeviationWithBiasCorrection(),
			}
			momentArr[i] = append(momentArr[i], moments)
		}
	}

	for i := 0; i < repeats; i++ {
		results := <-ch
		dfile.WriteString(fmt.Sprintf("%g,%g\n", results.ks, results.vd))
		for j := 0; j < maxl; j++ {
			momentArr[0][j].Increment(results.scovs[j])
			momentArr[1][j].Increment(results.rcovs[j])
			momentArr[2][j].Increment(results.xyPL[j])
			momentArr[3][j].Increment(results.xsysPL[j])
			momentArr[4][j].Increment(results.smXYPL[j])
		}

		if (i+1)%(repeats/100) == 0 {
			cfile, err := os.Create(fmt.Sprintf("%s_covs.csv", prefix))
			if err != nil {
				panic(err)
			}

			cfile.WriteString(fmt.Sprintf("#size: %d\n", size))
			cfile.WriteString(fmt.Sprintf("#length: %d\n", length))
			cfile.WriteString(fmt.Sprintf("#fragment: %d\n", fragment))
			cfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
			cfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
			cfile.WriteString(fmt.Sprintf("#sample: %d\n", sample))
			cfile.WriteString(fmt.Sprintf("#replicates: %d\n", i+1))
			cfile.WriteString("#dist, scov, rcov, xy, xsys, smxy_sd, scov_sd, rcov_sd, xy_sd, xsys_sd, smxy_sd\n")

			cfile.Close()
		}

	}
}
