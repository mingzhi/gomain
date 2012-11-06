package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"log"
	"os"
	"path"
	"strconv"
)

var (
	dir string // folder of files
	num int    // number of files
)

type Moments struct {
	Mean *desc.Mean
	Sd   *desc.StandardDeviation
}

func (m *Moments) Increment(v float64) {
	m.Sd.Increment(v)
	m.Mean.Increment(v)
}

func init() {
	flag.StringVar(&dir, "dir", "out", "folder of files")
	flag.IntVar(&num, "num", 1000, "number of files")

	flag.Parse()
}

func main() {
	momentArr := make([][]Moments, 5)
	for i := 0; i < num; i++ {
		filename := fmt.Sprintf("%s/covs_%d.csv", dir, i)

		f, err := os.Open(filename)
		if err != nil {
			log.Panic(err)
		}
		defer f.Close()

		reader := csv.NewReader(f)
		reader.Comment = '#'

		records, err := reader.ReadAll()
		if err != nil {
			log.Panic(err)
		}

		for j, rec := range records {
			if j >= len(momentArr[0]) {
				for k := 0; k < len(momentArr); k++ {
					m := Moments{
						Mean: desc.NewMean(),
						Sd:   desc.NewStandardDeviationWithBiasCorrection(),
					}
					momentArr[k] = append(momentArr[k], m)
				}
			}
			for k := 0; k < len(momentArr); k++ {
				v, err := strconv.ParseFloat(rec[k+1], 64)
				if err != nil {
					log.Panic(err)
				}
				momentArr[k][j].Increment(v)
			}
		}
	}

	filename := fmt.Sprintf("%s_covs.csv", path.Base(dir))
	f, err := os.Create(filename)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	for i := 0; i < len(momentArr[0]); i++ {
		f.WriteString(fmt.Sprintf("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
			i,
			momentArr[0][i].Mean.GetResult(),
			momentArr[1][i].Mean.GetResult(),
			momentArr[2][i].Mean.GetResult(),
			momentArr[3][i].Mean.GetResult(),
			momentArr[4][i].Mean.GetResult(),
			momentArr[0][i].Sd.GetResult(),
			momentArr[1][i].Sd.GetResult(),
			momentArr[2][i].Sd.GetResult(),
			momentArr[3][i].Sd.GetResult(),
			momentArr[4][i].Sd.GetResult(),
		))
	}
}
