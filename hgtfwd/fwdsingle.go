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
)

type Moment struct {
	Mean *desc.Mean
	Sd   *desc.StandardDeviation
}

func (m *Moment) Increment(d float64) {
	m.Mean.Increment(d)
	m.Sd.Increment(d)
}

var (
	size     int     // population size
	lens     int     // genome length
	mutation float64 // mutation rate per site per generation
	transfer float64 // transfer rate per site per generation
	frag     int     // fragment length to transfer
	maxl     int     // max distance to calculate
	gens     int     // total generations to sample after 10000 generation
	prefix   string  // prefix
	eqvGens  int     // generations to reach equiliquim
)

func init() {
	// register flags
	flag.IntVar(&size, "size", 1000, "population size")
	flag.IntVar(&lens, "genome", 1000, "genome length")
	flag.IntVar(&frag, "frag", 100, "fragment length to transfer")
	flag.IntVar(&maxl, "maxl", 200, "max distance to calculate")
	flag.IntVar(&gens, "gens", 10000, "total generations to sample after 10000 generation")
	flag.Float64Var(&transfer, "transfer", 1e-4, "transfer rate per site per generation")
	flag.Float64Var(&mutation, "mutation", 1e-4, "mutation rate per site per generation")
	flag.StringVar(&prefix, "prefix", "test", "prefix")
	flag.IntVar(&eqvGens, "eqv", 100000, "generations to reach equiliquim")

	// parse flags
	flag.Parse()
}

func main() {
	// create d file
	dfile, err := os.Create(fmt.Sprintf("%s_d.csv", prefix))
	if err != nil {
		log.Panic(err)
	}
	defer dfile.Close()

	// create population for simulation
	sp := fwd.NewSeqPop(size, lens, mutation, transfer, frag)

	// do 10000 generations for reaching equilibrium
	for i := 0; i < eqvGens; i++ {
		sp.Evolve()
	}

	// create moments
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

	// do sample generations
	for i := 0; i < gens; i++ {
		sp.Evolve()

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

		for j := 0; j < maxl; j++ {
			moments[0][j].Increment(scovs[j])
			moments[1][j].Increment(rcovs[j])
			moments[2][j].Increment(xyPL[j])
			moments[3][j].Increment(xsysPL[j])
			moments[4][j].Increment(smXYPL[j])
		}

		if (i+1)%(gens/100) == 0 {
			cfile, err := os.Create(fmt.Sprintf("%s_covs.csv", prefix))
			if err != nil {
				log.Panic(err)
			}

			cfile.WriteString(fmt.Sprintf("#size: %d\n", size))
			cfile.WriteString(fmt.Sprintf("#length: %d\n", lens))
			cfile.WriteString(fmt.Sprintf("#fragment: %d\n", frag))
			cfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
			cfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
			cfile.WriteString(fmt.Sprintf("#generations: %d\n", i+1))
			cfile.WriteString("dist, scov, rcov, xy, xsys, smxy_sd, scov_sd, rcov_sd, xy_sd, xsys_sd, smxy_sd\n")

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
			log.Printf("Finish %%%d\n", (i+1)/(gens/100))
		}

		dfile.WriteString(fmt.Sprintf("%g,%g\n", ks, vd))
	}
}
