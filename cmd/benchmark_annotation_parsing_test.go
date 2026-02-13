package cmd

import (
	"flag"
	"os"
	"strings"
	"testing"
)

var annofilePath = flag.String("path", "", "path to the input file")
var region = flag.String("region", "", "region of the chromosome to pull annotations for")

// This has to be called TestMain exactly to work
func TestMain(m *testing.M) {
	flag.Parse()

	os.Exit(m.Run())
}

func BenchmarkAnnoParser(b *testing.B) {

	keep_cols := "Consequence,Amino_acids,gnomADe_NFE_AF,gnomADg_NFE_AF,CLIN_SIG,PUBMED"
	keep_col_list := strings.Split(keep_cols, ",")

	parsed_region, _ := parse_region(*region)

	b.Logf("Running benchmarks")

	for b.Loop() {
		read_annotations(*annofilePath, keep_col_list, parsed_region)
	}
}
