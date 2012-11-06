package main

import (
	"flag"
	"fmt"
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
	ch := make(chan bool, ncpu)

	t0 := time.Now()

	for i := 0; i < ncpu; i++ {
		begin := i * repeats / ncpu
		end := (i + 1) * repeats / ncpu
		go simusome(begin, end, size, sample, length, mutation, transfer, fragment, prefix, ch)
	}

	for i := 0; i < ncpu; i++ {
		<-ch
	}

	t1 := time.Now()
	fmt.Println(t1.Sub(t0))
}

func simusome(begin, end, size, sample, length int, mutation, transfer float64, fragment int, prefix string, ch chan bool) {
	for i := begin; i < end; i++ {
		w := coals.NewWFPopulation(size, sample, length, mutation, transfer, fragment)
		simu(w, i, prefix)
	}
	ch <- true
}

func simu(w *coals.WFPopulation, idx int, prefix string) {
	dfile, err := os.Create(fmt.Sprintf("%sd_%d.csv", prefix, idx))
	if err != nil {
		panic(err)
	}
	defer dfile.Close()

	dfile.WriteString(fmt.Sprintf("#size: %d\n", w.Size))
	dfile.WriteString(fmt.Sprintf("#length: %d\n", w.GenomeLength))
	dfile.WriteString(fmt.Sprintf("#fragment: %d\n", w.TransferLength))
	dfile.WriteString(fmt.Sprintf("#mutation: %g\n", w.MutationRate))
	dfile.WriteString(fmt.Sprintf("#transfer: %g\n", w.TransferRate))
	dfile.WriteString(fmt.Sprintf("#sample: %d\n", w.SampleSize))
	dfile.WriteString("#ks, vd\n")

	w.Seed(idx)
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
	dfile.WriteString(fmt.Sprintf("%g,%g\n", ks, vd))

	scovs, rcovs, xyPL, xsysPL, smXYPL := cmatrix.CovCircle(maxl)

	covfile, err := os.Create(fmt.Sprintf("%scovs_%d.csv", prefix, idx))
	if err != nil {
		fmt.Println(err)
	}
	defer covfile.Close()

	covfile.WriteString(fmt.Sprintf("#size: %d\n", w.Size))
	covfile.WriteString(fmt.Sprintf("#length: %d\n", w.GenomeLength))
	covfile.WriteString(fmt.Sprintf("#fragment: %d\n", w.TransferLength))
	covfile.WriteString(fmt.Sprintf("#mutation: %g\n", w.MutationRate))
	covfile.WriteString(fmt.Sprintf("#transfer: %g\n", w.TransferRate))
	covfile.WriteString(fmt.Sprintf("#sample: %d\n", w.SampleSize))
	covfile.WriteString("#dist, scov, rcov, xy, xsys, smxy\n")

	for i := 0; i < maxl; i++ {
		covfile.WriteString(fmt.Sprintf("%d,%g,%g,%g,%g,%g\n",
			i,
			scovs[i],
			rcovs[i],
			xyPL[i],
			xsysPL[i],
			smXYPL[i],
		))
	}
}
