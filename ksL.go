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
	size := 1000       // population size = 1000
	length := 100000   // genome length = 100000
	mutation := 0.0001 // mutation rate = 1e-4
	transfer := 0.0    // transfer rate = 0
	fragment := 0      // transferred length = 0

	// construct a population
	pop := fwd.NewSeqPop(size, length, mutation, transfer, fragment)
	// simulation parameters
	samplesize := 100
	numofgen := 100

	// create files
	// file for length 100
	file1, err := os.Create("D_length_100.csv")
	if err != nil {
		panic(err)
	}
	defer file1.Close()
	// file for length 1000
	file2, err := os.Create("D_length_1000.csv")
	if err != nil {
		panic(err)
	}
	defer file2.Close()
	// file for length 10000
	file3, err := os.Create("D_length_10000.csv")
	if err != nil {
		panic(err)
	}
	defer file3.Close()
	// file for length 100000
	file4, err := os.Create("D_length_100000.csv")
	if err != nil {
		panic(err)
	}
	defer file4.Close()

	// do the simulation
	for i := 0; i < numofgen; i++ {
		pop.Evolve()

		// get #{samplesize} samples
		rSeq := rand.Perm(pop.Size) // get a permutation sequence
		samples := make([]fwd.Sequence, samplesize)
		for j := 0; j < samplesize; j++ {
			samples[j] = pop.Genomes[rSeq[j]]
		}

		// calculate the distance matrix for different lengthes
		dmatrix1 := [][]int{} // for length 100
		dmatrix2 := [][]int{} // for length 1000
		dmatrix3 := [][]int{} // for length 10000
		dmatrix4 := [][]int{} // for length 100000
		for j := 0; j < samplesize; j++ {
			for k := j + 1; k < samplesize; k++ {
				a := samples[j]
				b := samples[k]
				ds1 := []int{} // for length 100
				ds2 := []int{} // for length 1000
				ds3 := []int{} // for length 10000
				ds4 := []int{} // for length 100000
				for l := 0; l < len(a); l++ {
					if a[l] != b[l] {
						if l < 100 { // for length 100
							ds1 = append(ds1, l)
						}
						if l < 1000 { // for length 1000
							ds2 = append(ds2, l)
						}
						if l < 10000 { // for length 10000
							ds3 = append(ds3, l)
						}
						ds4 = append(ds4, l)
					}

				}
				dmatrix1 = append(dmatrix1, ds1)
				dmatrix2 = append(dmatrix2, ds2)
				dmatrix3 = append(dmatrix3, ds3)
				dmatrix4 = append(dmatrix4, ds4)
			}
		}
		// create cmatrix
		cmatrix1 := covs.NewCMatrix(samplesize, 100, dmatrix1)
		cmatrix2 := covs.NewCMatrix(samplesize, 1000, dmatrix2)
		cmatrix3 := covs.NewCMatrix(samplesize, 10000, dmatrix3)
		cmatrix4 := covs.NewCMatrix(samplesize, 100000, dmatrix4)
		// calculate ks and vard
		ks1, vd1 := cmatrix1.D()
		ks2, vd2 := cmatrix2.D()
		ks3, vd3 := cmatrix3.D()
		ks4, vd4 := cmatrix4.D()
		// write to file
		file1.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ks1, vd1))
		file2.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ks2, vd2))
		file3.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ks3, vd3))
		file4.WriteString(fmt.Sprintf("%d,%g,%g\n", pop.NumOfGen, ks4, vd4))
	}
	t1 := time.Now()
	log.Printf("End at: %v\n", t1)
	log.Printf("Duration: %v\n", t1.Sub(t0))
}
