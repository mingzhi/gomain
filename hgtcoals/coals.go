package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/hgt/coals"
	"github.com/mingzhi/hgt/covs"
	"math"
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
	means := make([][]*desc.Mean, 5)
	sds := make([][]*desc.StandardDeviation, 5)
	for i := 0; i < 5; i++ {
		means[i] = make([]*desc.Mean, maxl)
		sds[i] = make([]*desc.StandardDeviation, maxl)
		for j := 0; j < maxl; j++ {
			means[i][j] = desc.NewMean()
			sds[i][j] = desc.NewStandardDeviationWithBiasCorrection()
		}
	}

	dfile, err := os.Create(prefix + "d.csv")
	if err != nil {
		panic(err)
	}
	defer dfile.Close()

	dfile.WriteString(fmt.Sprintf("#size: %d\n", size))
	dfile.WriteString(fmt.Sprintf("#sample: %d\n", sample))
	dfile.WriteString(fmt.Sprintf("#length: %d\n", length))
	dfile.WriteString(fmt.Sprintf("#fragment: %d\n", fragment))
	dfile.WriteString(fmt.Sprintf("#repeats: %d\n", repeats))
	dfile.WriteString(fmt.Sprintf("#maxl: %d\n", maxl))
	dfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
	dfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
	dfile.WriteString("#ks, vd\n")

	runtime.GOMAXPROCS(runtime.NumCPU())

	t0 := time.Now()
	for c := 0; c < repeats; c++ {
		w := coals.NewWFPopulation(size, sample, length, mutation, transfer, fragment)
		w.Seed(c)
		w.Backtrace()
		seqs := w.Fortrace()
		diffmatrix := [][]int{}
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
		dfile.WriteString(fmt.Sprintf("%d,%g,%g\n", c+1, ks, vd))

		scovs, rcovs, xyPL, xsysPL, smXYPL := cmatrix.CovCircle(maxl)
		for l := 0; l < maxl; l++ {
			means[0][l].Increment(scovs[l])
			means[1][l].Increment(rcovs[l])
			means[2][l].Increment(xyPL[l])
			means[3][l].Increment(xsysPL[l])
			means[4][l].Increment(smXYPL[l])

			sds[0][l].Increment(scovs[l])
			sds[1][l].Increment(rcovs[l])
			sds[2][l].Increment(xyPL[l])
			sds[3][l].Increment(xsysPL[l])
			sds[4][l].Increment(smXYPL[l])
		}

		if (c+1)%(repeats/100) == 0 {
			t1 := time.Now()
			fmt.Printf("%d%%,%v\n", (c+1)/(repeats/100), t1.Sub(t0))
		}
	}

	covfile, err := os.Create(prefix + "covs.csv")
	if err != nil {
		fmt.Println(err)
	}
	defer covfile.Close()

	covfile.WriteString(fmt.Sprintf("#size: %d\n", size))
	covfile.WriteString(fmt.Sprintf("#sample: %d\n", sample))
	covfile.WriteString(fmt.Sprintf("#length: %d\n", length))
	covfile.WriteString(fmt.Sprintf("#fragment: %d\n", fragment))
	covfile.WriteString(fmt.Sprintf("#repeats: %d\n", repeats))
	covfile.WriteString(fmt.Sprintf("#maxl: %d\n", maxl))
	covfile.WriteString(fmt.Sprintf("#mutation: %g\n", mutation))
	covfile.WriteString(fmt.Sprintf("#transfer: %g\n", transfer))
	covfile.WriteString("#dist, scov, rcov, xy, xsys, smxy, scov_sd, rcov_sd, xy_sd, xsys_sd, smxy_sd\n")

	for i := 0; i < maxl; i++ {
		covfile.WriteString(fmt.Sprintf("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
			i,
			means[0][i].GetResult(),
			means[1][i].GetResult(),
			means[2][i].GetResult(),
			means[3][i].GetResult(),
			means[4][i].GetResult(),
			sds[0][i].GetResult()/math.Sqrt(float64(repeats)),
			sds[1][i].GetResult()/math.Sqrt(float64(repeats)),
			sds[2][i].GetResult()/math.Sqrt(float64(repeats)),
			sds[3][i].GetResult()/math.Sqrt(float64(repeats)),
			sds[4][i].GetResult()/math.Sqrt(float64(repeats)),
		))
	}
}
