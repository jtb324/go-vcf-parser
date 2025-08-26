package files

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

func read_header(header_line string) map[int]string {
	mapping_dict := make(map[int]string)

	split_header := strings.Split(strings.TrimSpace(header_line), "\t")

	for indx, colname := range split_header {
		mapping_dict[indx] = colname
	}

	return mapping_dict
}

type VcfStreamer struct {
	Scanner             *bufio.Scanner
	Header_col_mappings map[int]string
	Next_line           string
}

func NewVcfStreamer(bufferSize int) *VcfStreamer {
	buf := make([]byte, 0, bufferSize)
	scanner := bufio.NewScanner(os.Stdin)
	scanner.Buffer(buf, bufferSize)
	return &VcfStreamer{Scanner: scanner}
}

// This stream will have info lines, the header lines, and then it will get to the variant.
// We are going to process the info lines and the header lines here
func (vcfStreamer *VcfStreamer) Initialize() error {
	for vcfStreamer.Scanner.Scan() {
		line := vcfStreamer.Scanner.Text()
		if strings.Contains(line, "##") {
			continue
		} else if strings.Contains(line, "#CHROM") {
			vcfStreamer.Header_col_mappings = read_header(line)
		} else {
			vcfStreamer.Next_line = line
			break
		}
	}
	if vcfStreamer.Scanner.Err() != nil {
		return fmt.Errorf("the following error was encountered while trying to read through the vcf info lines and header lines: %s", vcfStreamer.Scanner.Err())
	}
	return nil
}

func (vcfStreamer *VcfStreamer) ReadNextLine() {
	if vcfStreamer.Scanner.Scan() {
		vcfStreamer.Next_line = vcfStreamer.Scanner.Text()
	} else {
		vcfStreamer.Next_line = ""
	}
}

func (vcfStreamer *VcfStreamer) CheckErrs() error {
	if vcfStreamer.Scanner.Err() != nil {
		return fmt.Errorf("encountered the following error while attempting to parse the input vcf file stream: %s", vcfStreamer.Scanner.Err())
	}
	return nil
}
