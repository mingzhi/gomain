// in this script, I want to determine the impact of genome length.
package main

import (
	"fmt"
	covs "github.com/mingzhi/hgt/covs"
	fwd "github.com/mingzhi/hgt/fwd3"
	"log"
	"math/rand"
	"os"
	"time"
)

func main() {
	t0 := time.Now() // start time
	log.Printf("Start up at: %v\n", t0)

	// population paramters
	length := 1000     // genome length = 100000
	mutation := 0.0001 // mutation rate = 1e-4
	transfer := 0.0    // transfer rate = 0
	fragment := 0      // transferred length = 0

	// simulation parameters
	samplesize := 100
	numofgen := 100000

	sizearray := []int{100, 1000, 10000, 100000}
	for _, size := range sizearray {
		file, err := os.Create(fmt.Sprintf("D_size_%d.csv", size))
		if err != nil {
			panic(err)
		}
		pop := fwd.NewSeqPop(size, length, mutation, transfer, fragment)
		for i := 0; i < numofgen; i++ {
			pop.Evolve()
			// get #{samplesize} samples
			rSeq := rand.Perm(pop.Size) // get a permutation sequence
			samples := make([]fwd.Sequence, samplesize)
			for j := 0; j < samplesize; j++ {
				samples[j] = pop.Genomes[rSeq[j]]
			}

			// calculate the distance matrix for different lengthes
			dmatrix := [][]int{}
			for j := 0; j < samplesize; j++ {
				for k := j + 1; k < samplesize; k++ {
					a := samples[j]
					b := samples[k]
					ds := []int{}
					for l := 0; l < len(a); l++ {
						if a[l] != b[l] {
							ds = append(ds, l)
						}

					}
					dmatrix = append(dmatrix, ds)
				}
			}
			// create cmatrix
			cmatrix := covs.NewCMatrix(samplesize, pop.Length, dmatrix)
			// calculate ks and vard
			ks, vd := cmatrix.D()
			// write to file
			file.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ks, vd))

			// printing the process
			if pop.NumOfGen%1000 == 0 {
				t2 := time.Now()
				fmt.Printf("Number of Generation with size = %d: %d, time used = %v\n", pop.Size, pop.NumOfGen, t2.Sub(t0))
			}
		}

		file.Close()
	}

	t1 := time.Now()
	log.Printf("End at: %v\n", t1)
	log.Printf("Duration: %v\n", t1.Sub(t0))
}
