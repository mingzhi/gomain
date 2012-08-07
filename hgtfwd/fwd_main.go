package main

import (
	"bitbucket.org/mingzhi/gsl/randist"
	"bitbucket.org/mingzhi/hgt/covs"
	"bitbucket.org/mingzhi/hgt/fwd"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
)

var (
	size         int
	length       int
	fragment     int
	mutationRate float64
	transferRate float64
	time         int
	sampleSize   int
	stepSize     int
	maxL         int
	prex         string
	dir          string
	pop          *fwd.SeqPopulation
)

func init() {
	//assign flags
	flag.IntVar(&size, "size", 1000, "usage: -size 1000")
	flag.IntVar(&length, "length", 100000, "usage: -length 1000")
	flag.IntVar(&time, "time", 1000, "usage: -time 1000")
	flag.IntVar(&sampleSize, "samplesize", 100, "usage: -samplesize 1000")
	flag.IntVar(&stepSize, "stepsize", 1000, "usage: -stepsize 1000")
	flag.IntVar(&maxL, "maxl", 1000, "usage: -maxl 1000")
	flag.IntVar(&fragment, "fragment", 1000, "usage: -fragment 1000")
	flag.Float64Var(&mutationRate, "mutation", 1e-8, "usage: -mutation 1e-8")
	flag.Float64Var(&transferRate, "transfer", 0, "usage: -transfer 1e-9")
	flag.StringVar(&prex, "prex", "default", "usage: -prex default")
	flag.StringVar(&dir, "out", "out", "usage: -dir out")

	// parse flags
	flag.Parse()
	log.SetPrefix(prex + ":")

	// set random number generator
	rng := randist.NewRNG(randist.MT19937)

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

	stepNum := time / stepSize
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

		// calculate covs
		covfilestr := dir + "/" + prex + fmt.Sprintf("_covs_%d.txt", pop.Time)
		covfile, err := os.Create(covfilestr)
		if err != nil {
			log.Fatalf("Can not create file: %s, %v", covfilestr, err)
		}
		scovs, rcovs, xyPL, xsysPL, smXYPL := c.CalculateCovs(maxL)
		for j := 0; j < maxL; j++ {
			covfile.WriteString(fmt.Sprintf("%d\t%g\t%g\t%g\t%g\t%g\n", j, scovs[j], rcovs[j], xyPL[j], xsysPL[j], smXYPL[j]))
		}
		covfile.Close()
	}
}
