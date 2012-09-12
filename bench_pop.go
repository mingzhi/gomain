package main

import (
	"fmt"
	fwd "github.com/mingzhi/hgt/fwd1"
	"time"
)

func main() {
	t0 := time.Now()

	size := 1000
	length := 1000
	mutation := 0.001
	transfer := 0.0
	fragment := 0

	pop := fwd.NewSeqPop(size, length, mutation, transfer, fragment)

	n := 1000
	for i := 0; i < n; i++ {
		pop.Evolve()
	}

	t1 := time.Now()
	t := t1.Sub(t0)
	fmt.Println(t)
}
