package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"

	gzip "github.com/klauspost/compress/gzip"

	_ "net/http/pprof"
)

var anno_file = "/Users/jtbaker/code_projects/go-phers-parser/test_data/agd250k_chr22_vep_annotations.vcf.gz"

func parse_annotations_no_gzip(annotation_filepath string, buffer_size int) {
	// we need to open the file, pass to a gzip reader and then pass to a bufio Scanner
	fmt.Println("Parsing annotations")
	anno_fh, anno_err := os.Open(annotation_filepath)

	if anno_err != nil {
		fmt.Errorf("encountered the following error while attempting to open the annotation file:\n %s", anno_err)
	}

	defer anno_fh.Close()

	buf := make([]byte, 5120*5120)

	anno_scanner := bufio.NewScanner(anno_fh)

	anno_scanner.Buffer(buf, 5120*5120)

	// create a byte reader

	for anno_scanner.Scan() {
		reader := bytes.NewBuffer(anno_scanner.Bytes())

		gzreader, _ := gzip.NewReader(reader)

		if val, err := io.ReadAll(gzreader); err == nil {
			fmt.Println(string(val))
		} else {
			fmt.Printf("%s\n", err)
		}
		// if output, e2 := io.ReadAll(gzreader); e2 == nil {
		// 	fmt.Println(string(output))
		// 	break
		// } else {
		// 	fmt.Printf("Error in the io.ReadAll: %s\n", e2)
		// 	break
		// }
		// if val, err := gzreader.Read(b); err == nil {
		// 	fmt.Println(val)
		// 	fmt.Println(b)
		// 	break
		// } else {
		// 	fmt.Printf("%s\n", err)
		// }

		break
	}
	if anno_scanner.Err() != nil {
		fmt.Errorf("There was an error while reading in the annotation_filepath: %s\n", anno_scanner.Err())
	}
}

func main() {
	// go func() {
	// 	log.Println(http.ListenAndServe("localhost:6060", nil))
	// }()
	parse_annotations_no_gzip(anno_file, 4096)
}
