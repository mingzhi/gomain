package main

import (
	"fmt"
	covs "github.com/mingzhi/hgt/covs"
	fwd "github.com/mingzhi/hgt/fwd1"
	"os"
	"runtime"
)

func main() {
	// use all cpus
	runtime.GOMAXPROCS(runtime.NumCPU())

	// population paramters
	length := 1000      // genome length = 1000
	mutation := 0.00001 // mutation rate = 1e-5
	transfer := 0.0     // transfer rate = 0
	fragment := 0       // transferred length = 0

	// simulation parameters
	samplesize := 100
	numofgen := 100000

	// population size array
	sizearray := []int{100, 1000, 10000}
	for _, size := range sizearray {
		// population
		pop := fwd.NewSeqPop(size, length, mutation, transfer, fragment)

		// result files
		ksname := fmt.Sprintf("ks_size_%d.csv", size)
		ksfile, err := os.Create(ksname)
		if err != nil {
			panic(err)
		}
		vdname := fmt.Sprintf("vd_size_%d.csv", size)
		vdfile, err := os.Create(vdname)
		if err != nil {
			panic(err)
		}

		// info file
		infoname := fmt.Sprintf("ks_size_%d.txt", size)
		infofile, err := os.Create(infoname)
		if err != nil {
			panic(err)
		}
		infofile.WriteString(fmt.Sprintf("Size = %d\n", pop.Size))
		infofile.WriteString(fmt.Sprintf("Length = %d\n", pop.Length))
		infofile.WriteString(fmt.Sprintf("Mutation = %g\n", pop.Mutation))

		// population evolve
		for i := 0; i < numofgen; i++ {
			pop.Evolve()
			sample := fwd.Sample(pop.Genomes, samplesize)
			dmatrix := fwd.GenerateDistanceMatrix(sample)
			cmatrix := covs.NewCMatrix(samplesize, length, dmatrix)
			ks, vd := cmatrix.D()
			ksfile.WriteString(fmt.Sprintf("%g\n", ks))
			vdfile.WriteString(fmt.Sprintf("%g\n", vd))
		}

		// close files
		ksfile.Close()
		vdfile.Close()
		infofile.Close()
	}
}
