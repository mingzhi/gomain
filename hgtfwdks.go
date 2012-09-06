package main

import (
	"fmt"
	"github.com/mingzhi/gomath/random"
	"github.com/mingzhi/hgt/covs"
	"github.com/mingzhi/hgt/fwd1"
	"log"
	"math/rand"
	"os"
	"runtime"
)

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())

	length := 1000
	mutation := 1e-5
	transfer := 0.0
	fragment := 0
	src := random.NewLockedSource(rand.NewSource(1))
	time := 100
	step := 10000
	sampleSize := 100

	sizearray := []int{100, 1000, 10000, 100000}
	for _, size := range sizearray {
		ksfilestr := fmt.Sprintf("ks_%d.csv", size)
		ksfile, err := os.Create(ksfilestr)
		if err != nil {
			log.Fatalf("Can not create file: %s, %v", ksfilestr, err)
		}
		pop := fwd1.NewSeqPop(size, length, mutation, transfer, fragment, src)
		for i := 0; i < step; i++ {
			for t := 0; t < time; t++ {
				pop.Evolve()
			}
			samples := fwd1.Sample(pop.Genomes, sampleSize)
			dmatrix := fwd1.GenerateDistanceMatrix(samples)
			cmatrix := covs.NewCMatrix(sampleSize, length, dmatrix)
			ks := cmatrix.KS()
			ksfile.WriteString(fmt.Sprintf("%d,%.10f\n", pop.NumOfGen, ks))
		}
		ksfile.Close()
	}
}
