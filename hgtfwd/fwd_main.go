package main

import (
	"bitbucket.org/mingzhi/chart/render"
	"bitbucket.org/mingzhi/gomath/stat/desc"
	"bitbucket.org/mingzhi/gsl/randist"
	"bitbucket.org/mingzhi/hgt/covs"
	"bitbucket.org/mingzhi/hgt/fwd"
	"flag"
	"fmt"
	"github.com/vdobler/chart"
	"log"
	"os"
	"runtime"
	"time"
)

var (
	size         int
	length       int
	fragment     int
	mutationRate float64
	transferRate float64
	etime        int
	sampleSize   int
	stepSize     int
	maxL         int
	prex         string
	dir          string
	pop          *fwd.SeqPopulation
	tnow         time.Time
)

func init() {
	//assign flags
	flag.IntVar(&size, "size", 1000, "usage: -size 1000")
	flag.IntVar(&length, "length", 100000, "usage: -length 1000")
	flag.IntVar(&etime, "time", 1000, "usage: -time 1000")
	flag.IntVar(&sampleSize, "samplesize", 100, "usage: -samplesize 1000")
	flag.IntVar(&stepSize, "stepsize", 1000, "usage: -stepsize 1000")
	flag.IntVar(&maxL, "maxl", 1000, "usage: -maxl 1000")
	flag.IntVar(&fragment, "fragment", 1000, "usage: -fragment 1000")
	flag.Float64Var(&mutationRate, "mutation", 1, "usage: -mutation 1e-8")
	flag.Float64Var(&transferRate, "transfer", 0, "usage: -transfer 1e-9")
	flag.StringVar(&prex, "prex", "default", "usage: -prex default")
	flag.StringVar(&dir, "out", "out", "usage: -dir out")

	// parse flags
	flag.Parse()
	log.SetPrefix(prex + ":")

	// get start time
	tnow = time.Now()
	log.Printf("Begin at %v\n", tnow.Format("Mon Jan 2 15:04:05"))

	// set random number generator
	rng := randist.NewRNG(randist.MT19937)
	seed := tnow.Nanosecond()
	rng.SetSeed(seed)
	log.Println("Set random seed:", seed)

	// init population
	pop = fwd.NewSeqPopulation(size, length, mutationRate, transferRate, fragment, rng)
	//log.Println("Population: ", pop)
	log.Println("Population initialized.")

	// determine how many cpus that we can use
	ncpu := runtime.NumCPU()
	runtime.GOMAXPROCS(ncpu)
	log.Println("Number of CPU used: ", ncpu)
}

func main() {
	// ks file
	ksfilestr := dir + "/" + prex + "_ks.txt"
	ksfile, err := os.Create(ksfilestr)
	if err != nil {
		log.Fatalf("Can not create file: %s, %v", ksfilestr, err)
	}
	defer ksfile.Close()

	ksMean := desc.NewMean()
	ksVar := desc.NewVarianceWithBiasCorrection()

	vardMean := desc.NewMean()
	vardVar := desc.NewVarianceWithBiasCorrection()

	scovsMeanArray := make([]*desc.Mean, maxL)
	scovsVarArray := make([]*desc.Variance, maxL)

	rcovsMeanArray := make([]*desc.Mean, maxL)
	rcovsVarArray := make([]*desc.Variance, maxL)

	xyPLMeanArray := make([]*desc.Mean, maxL)
	xyPLVarArray := make([]*desc.Variance, maxL)

	xsysPLMeanArray := make([]*desc.Mean, maxL)
	xsysPLVarArray := make([]*desc.Variance, maxL)

	smXYPLMeanArray := make([]*desc.Mean, maxL)
	smXYPLVarArray := make([]*desc.Variance, maxL)

	for i := 0; i < maxL; i++ {
		scovsMeanArray[i] = desc.NewMean()
		scovsVarArray[i] = desc.NewVarianceWithBiasCorrection()

		rcovsMeanArray[i] = desc.NewMean()
		rcovsVarArray[i] = desc.NewVarianceWithBiasCorrection()

		xyPLMeanArray[i] = desc.NewMean()
		xyPLVarArray[i] = desc.NewVarianceWithBiasCorrection()

		xsysPLMeanArray[i] = desc.NewMean()
		xsysPLVarArray[i] = desc.NewVarianceWithBiasCorrection()

		smXYPLMeanArray[i] = desc.NewMean()
		smXYPLVarArray[i] = desc.NewVarianceWithBiasCorrection()

	}

	stepNum := etime / stepSize
	for i := 0; i < stepNum; i++ {
		// simulate
		pop.Evolute(stepSize, fwd.K2PSubstitution)

		// create distance matrix
		dmatrix := pop.GenerateDistanceMatrix(sampleSize)
		c := covs.NewCMatrix(dmatrix, sampleSize, pop.Length)

		// calculate ks
		ks := c.CalculateKS()
		vard := c.CalculateVarD()
		ksfile.WriteString(fmt.Sprintf("%d\t%g\t%g\n", pop.Time, ks, vard))
		ksMean.Increment(ks)
		ksVar.Increment(ks)
		vardMean.Increment(vard)
		vardVar.Increment(vard)
		log.Printf("KsMean: %g, KsVar: %g, VardMean: %g, VardVar: %g\n", ksMean.GetResult(), ksVar.GetResult(), vardMean.GetResult(), vardVar.GetResult())

		// calculate covs
		scovs, rcovs, xyPL, xsysPL, smXYPL := c.CalculateCovs(maxL)
		for j := 0; j < maxL; j++ {
			scovsMeanArray[j].Increment(scovs[j])
			scovsVarArray[j].Increment(scovs[j])

			rcovsMeanArray[j].Increment(rcovs[j])
			rcovsVarArray[j].Increment(rcovs[j])

			xyPLMeanArray[j].Increment(xyPL[j])
			xyPLVarArray[j].Increment(xyPL[j])

			xsysPLMeanArray[j].Increment(xsysPL[j])
			xsysPLVarArray[j].Increment(xsysPL[j])

			smXYPLMeanArray[j].Increment(smXYPL[j])
			smXYPLVarArray[j].Increment(smXYPL[j])
		}

		svger := render.NewSVG(dir+"/"+prex+"_covs", 1, 5, 800, 200)
		scovMeans := make([]float64, maxL-1)
		rcovMeans := make([]float64, maxL-1)
		xyPLMeans := make([]float64, maxL-1)
		xsysPLMeans := make([]float64, maxL-1)
		smXYPLMeans := make([]float64, maxL-1)
		xs := make([]float64, maxL-1)
		for j := 1; j < maxL; j++ {
			xs[j-1] = float64(j)
			scovMeans[j-1] = scovsMeanArray[j].GetResult()
			rcovMeans[j-1] = rcovsMeanArray[j].GetResult()
			xyPLMeans[j-1] = xyPLMeanArray[j].GetResult()
			xsysPLMeans[j-1] = xsysPLMeanArray[j].GetResult()
			smXYPLMeans[j-1] = smXYPLMeanArray[j].GetResult()
		}
		pl := chart.ScatterChart{Title: "scovs"}
		pl.AddDataPair("struc covs", xs, scovMeans, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)

		pl = chart.ScatterChart{Title: "rcovs"}
		pl.AddDataPair("rate covs", xs, rcovMeans, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)

		pl = chart.ScatterChart{Title: "xyPL"}
		pl.AddDataPair("xyPL", xs, xyPLMeans, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)

		pl = chart.ScatterChart{Title: "xsysPL"}
		pl.AddDataPair("xsysPL", xs, xsysPLMeans, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)

		pl = chart.ScatterChart{Title: "smXYPL"}
		pl.AddDataPair("smXYPL", xs, smXYPLMeans, chart.PlotStyleLines, chart.Style{Symbol: '+', SymbolColor: "#0000ff", LineStyle: chart.SolidLine})
		svger.Plot(&pl)

		svger.Close()

		// write to file
		covsfilestr := dir + "/" + prex + "_covs.txt"
		covsfile, err := os.Create(covsfilestr)
		if err != nil {
			log.Fatalf("Can not create file: %s, %v", covsfilestr, err)
		}
		for j := 0; j < maxL; j++ {
			covsfile.WriteString(
				fmt.Sprintf("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					j,
					scovsMeanArray[j].GetResult(), scovsVarArray[j].GetResult(),
					rcovsMeanArray[j].GetResult(), rcovsVarArray[j].GetResult(),
					xyPLMeanArray[j].GetResult(), xyPLVarArray[j].GetResult(),
					xsysPLMeanArray[j].GetResult(), xsysPLVarArray[j].GetResult(),
					smXYPLMeanArray[j].GetResult(), smXYPLVarArray[j].GetResult(),
				),
			)
		}
		covsfile.Close()
	}

	tnow = time.Now()
	log.Printf("Done at %v\n", tnow)
}
