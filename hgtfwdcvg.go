package main

import (
	"flag"
	"github.com/mingzhi/chart/render"
	"github.com/mingzhi/gomath/random"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/covs"
	"github.com/mingzhi/hgt/fwd1"
	"github.com/vdobler/chart"
	"math/rand"
	"runtime"
)

var (
	size       int          // population size
	length     int          // genome length
	fragment   int          // transferred fragment
	mutation   float64      // mutation rate
	transfer   float64      // transfer rate
	ngen       int          // number of generation we want to evolve
	sampleSize int          // sample size
	pop        *fwd1.SeqPop // population
	fname      string       // file name
)

func init() {
	//assign flags
	flag.IntVar(&size, "size", 1000, "usage: -size 1000")
	flag.IntVar(&length, "length", 1000, "usage: -length 1000")
	flag.IntVar(&ngen, "ngen", 1000, "usage: -ngen 1000")
	flag.IntVar(&sampleSize, "samplesize", 100, "usage: -samplesize 1000")
	flag.IntVar(&fragment, "fragment", 100, "usage: -fragment 1000")
	flag.Float64Var(&mutation, "mutation", 1e-5, "usage: -mutation 1e-5")
	flag.Float64Var(&transfer, "transfer", 0.0, "usage: -transfer 1e-9")
	flag.StringVar(&fname, "fname", "ks", "usage: -fname ks")

	// parse flags
	flag.Parse()

	// create random source (locked)
	src := random.NewLockedSource(rand.NewSource(1))

	// construct a population
	pop = fwd1.NewSeqPop(size, length, mutation, transfer, fragment, src)
	// use all the available CPUs
	runtime.GOMAXPROCS(runtime.NumCPU())
}

// the converge of KS equalibrium.
func main() {
	ksarray := []float64{} // store ks
	ngarray := []float64{} // store generation number
	// do evolution
	for i := 0; i < ngen; i++ {
		pop.Evolve()
		// select 10 samples and averge
		mean := desc.NewMean()
		num := 10
		for j := 0; j < num; j++ {
			sample := fwd1.Sample(pop.Genomes, sampleSize)
			dmatrix := fwd1.GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(sampleSize, pop.Length, dmatrix)
			ks := cmatrix.KS()
			mean.Increment(ks)
		}

		ksarray = append(ksarray, mean.GetResult())
		ngarray = append(ngarray, float64(pop.NumOfGen))
	}

	// draw
	svger := render.NewSVG(fname, 1, 1, 800, 200)
	pl := chart.ScatterChart{Title: "KS"}
	pl.AddDataPair("KS", ngarray, ksarray, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
	svger.Plot(&pl)
	svger.Close()
}
