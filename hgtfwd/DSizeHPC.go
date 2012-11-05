// in this script, I want to determine the impact of genome length.
package main

import (
	"fmt"
	covs "github.com/mingzhi/hgt/covs"
	fwd "github.com/mingzhi/hgt/fwd3"
	"log"
	"math/rand"
	"os"
	"runtime"
	"time"
)

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())
	t0 := time.Now() // start time
	log.Printf("Start up at: %v\n", t0)

	// population paramters
	length := 1000     // genome length = 100000
	transfer := 0.0    // transfer rate = 0
	fragment := 0      // transferred length = 0
	mutation := 0.0001 // mutation rate = 1e-4

	// simulation parameters
	samplesize := 100
	numofgen := 100000

	sizearray := []int{100, 1000, 10000, 100000}
	ch := make(chan bool)
	for _, size := range sizearray {
		fname := fmt.Sprintf("D_size_%d.csv", size)
		go evolve(fname, size, length, mutation, transfer, fragment, samplesize, numofgen, ch)
	}

	for i := 0; i < len(sizearray); i++ {
		<-ch
	}

	t1 := time.Now()
	log.Printf("End at: %v\n", t1)
	log.Printf("Duration: %v\n", t1.Sub(t0))
}

func evolve(fname string, size, length int, mutation, transfer float64, fragment, samplesize, numofgen int, ch chan bool) {

	file, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer file.Close()
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
	}
	ch <- true
}
