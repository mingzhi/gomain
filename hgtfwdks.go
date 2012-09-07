package main

import (
	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"code.google.com/p/plotinum/vg"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/covs"
	fwd "github.com/mingzhi/hgt/fwd1"
	"os"
	"runtime"
)

var (
	size       int         // population size
	length     int         // genome length
	fragment   int         // transferred fragment
	mutation   float64     // mutation rate
	transfer   float64     // transfer rate
	ngen       int         // number of generation we want to evolve
	sampleSize int         // sample size
	sampleTime int         // sample times
	pop        *fwd.SeqPop // population
	fname      string      // file name
)

func init() {
	//assign flags
	flag.IntVar(&size, "size", 1000, "usage: -size 1000")
	flag.IntVar(&length, "length", 1000, "usage: -length 1000")
	flag.IntVar(&ngen, "ngen", 1000, "usage: -ngen 1000")
	flag.IntVar(&sampleSize, "samplesize", 100, "usage: -samplesize 1000")
	flag.IntVar(&sampleTime, "sampletime", 1, "usage: -sampletime 10")
	flag.IntVar(&fragment, "fragment", 100, "usage: -fragment 1000")
	flag.Float64Var(&mutation, "mutation", 1e-5, "usage: -mutation 1e-5")
	flag.Float64Var(&transfer, "transfer", 0.0, "usage: -transfer 1e-9")
	flag.StringVar(&fname, "fname", "ks", "usage: -fname ks")

	// parse flags
	flag.Parse()

	// construct a population
	pop = fwd.NewSeqPop(size, length, mutation, transfer, fragment)
	// use all the available CPUs
	runtime.GOMAXPROCS(runtime.NumCPU())
}

func main() {
	// create file storing ks and vard
	f, err := os.Create(fname + ".csv")
	if err != nil {
		panic(err)
	}
	defer f.Close()

	ksarray := []float64{} // store ks
	vdarray := []float64{} // store VarD
	ngarray := []float64{} // store generation number

	// do evolution
	for i := 0; i < ngen; i++ {
		pop.Evolve()
		// we make 10 samples and average
		ksmean := desc.NewMean()
		vdmean := desc.NewMean()
		for j := 0; j < sampleTime; j++ {
			sample := fwd.Sample(pop.Genomes, sampleSize)
			dmatrix := fwd.GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(sampleSize, pop.Length, dmatrix)
			ks, vard := cmatrix.D()
			ksmean.Increment(ks)
			vdmean.Increment(vard)
		}
		// write to csv
		f.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ksmean.GetResult(), vdmean.GetResult()))
		if (i+1)%1000 == 0 {
			fmt.Println("Generation: ", pop.NumOfGen, ksmean.GetResult(), vdmean.GetResult())
		}

		// store array for draw
		ksarray = append(ksarray, ksmean.GetResult())
		vdarray = append(vdarray, vdmean.GetResult())
		ngarray = append(ngarray, float64(pop.NumOfGen))
	}

	// draw
	/*
		svger := render.NewSVG(fname, 1, 2, 800, 200)
		pl := chart.ScatterChart{Title: "KS"}
		pl.AddDataPair("KS", ngarray, ksarray, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)
		pl = chart.ScatterChart{Title: "VarD"}
		pl.AddDataPair("VarD", ngarray, vdarray, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)
		svger.Close()
	*/
	// create a new plot, set its title and axis label
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "KS"

	// save population
	jf, err := os.Create(fname + ".json")
	if err != nil {
		panic(err)
	}
	defer jf.Close()
	b, err := json.Marshal(pop)
	if err != nil {
		panic(err)
	}
	jf.Write(b)
}
