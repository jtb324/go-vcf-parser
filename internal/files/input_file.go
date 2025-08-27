package files

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	gzip "github.com/klauspost/pgzip"
)

type Scanner interface {
	Scan() bool
	Text() string
	Err() error
}

type FileObj interface {
	Close() error
}

type FileReader struct {
	Filename        string
	FileScanner     *bufio.Scanner
	Err             error
	Header_col_indx map[string]int
	Header_Found    bool
	Col_count       int
	Handles         []io.Closer
}

func (fr FileReader) CheckErrors() {
	if errors.Is(fr.Err, os.ErrNotExist) {
		fmt.Printf("The file %s was not found.\n", fr.Filename)
	} else if errors.Is(fr.Err, os.ErrPermission) {
		fileStats, _ := os.Stat(fr.Filename)
		fmt.Printf("Couldn't open the file %s. Current file permissions are %s\n", fr.Filename, fileStats.Mode().Perm())
	} else {
		fmt.Printf("Encountered the following error while trying to open the file %s\n %s\n", fr.Filename, fr.Err)
	}
	os.Exit(1)
}

func mapHeader(header_line string) (map[string]int, int) {
	column_mappings := make(map[string]int)

	column_list := strings.Split(strings.TrimSpace(header_line), "\t")

	for indx, value := range column_list {
		column_mappings[value] = indx
	}

	return column_mappings, len(column_list)
}

func (fr *FileReader) ParseHeader(headerIdentified string) error {
	for fr.FileScanner.Scan() {
		line := fr.FileScanner.Text()
		if strings.Contains(line, headerIdentified) {
			col_indx, col_count := mapHeader(line)
			// We will need to use the column indices and the col count later
			fr.Header_col_indx = col_indx
			fr.Col_count = col_count
			// We also need to update that the header was found
			fr.Header_Found = true
			break
		}
	}
	if fr.FileScanner.Err() != nil {
		return fr.FileScanner.Err()
	}
	return nil
}

// Handle the creation of the file reader and the creation of a bufio.Scanner
func MakeCompressedFileReader(filename string, buffersize int) *FileReader {
	handles := make([]io.Closer, 2)

	fh, open_err := os.Open(filename)

	if open_err != nil {
		return &FileReader{Filename: filename, FileScanner: nil, Err: fmt.Errorf("encountered the following error while opening the file: %w", open_err), Handles: handles, Header_Found: false}
	}

	handles[0] = fh

	gh, gzip_err := gzip.NewReader(fh)

	if gzip_err != nil {
		return &FileReader{Filename: filename, FileScanner: nil, Err: fmt.Errorf("encountered the following error while trying to decompress the file: %w", gzip_err), Handles: handles, Header_Found: false}
	}

	handles[1] = gh

	buf := make([]byte, 0, buffersize)

	scanner := bufio.NewScanner(gh)

	scanner.Buffer(buf, buffersize)

	return &FileReader{Filename: filename, FileScanner: scanner, Err: nil, Handles: handles, Header_Found: false}
}

// Handle the creation of the file reader and the creation of a bufio.Scanner
func MakeFileReader(filename string, buffersize int) *FileReader {
	handles := make([]io.Closer, 1)
	fh, open_err := os.Open(filename)

	if open_err != nil {
		return &FileReader{Filename: filename, FileScanner: nil, Err: fmt.Errorf("encountered the following error while opening the file: %w", open_err), Handles: handles, Header_Found: false}
	}

	handles[0] = fh

	buf := make([]byte, 0, buffersize)

	scanner := bufio.NewScanner(fh)

	scanner.Buffer(buf, buffersize)

	return &FileReader{Filename: filename, FileScanner: scanner, Err: nil, Handles: handles, Header_Found: false}
}

func MakeStreamReader(buffersize int) *VCFReader {
	buf := make([]byte, 0, buffersize)

	stdin_streamer := bufio.NewScanner(os.Stdin)

	stdin_streamer.Buffer(buf, buffersize)

	fileReader := FileReader{
		Filename:    "standard input",
		FileScanner: stdin_streamer,
		Err:         nil,
		Handles:     nil,
	}

	return &VCFReader{FileReader: fileReader}
}

type VCFReader struct {
	FileReader
	SampleMapping    map[int]string
	SampleExclusions []string // Sometimes in VCF files there are samples that we want to ignore (reference panel samples or invalid samples). This attribute will help us ignore them
}

func (vcfReader *VCFReader) ParseHeader(header_identifier string) error {
	for vcfReader.FileScanner.Scan() {
		line := vcfReader.FileScanner.Text()
		if strings.Contains(line, header_identifier) {
			col_indx, col_count := mapHeader(line)
			// We will need to use the column indices and the col count later
			vcfReader.Header_col_indx = col_indx
			vcfReader.Col_count = col_count
			// Now we also have to map the sample ids where the key is the indx and the value is the column label
			vcfReader.SampleMapping = mapSamples(line, vcfReader.SampleExclusions)
			// We also need to update that the header was found
			vcfReader.Header_Found = true
			break
		}
	}
	if vcfReader.FileScanner.Err() != nil {
		return vcfReader.FileScanner.Err()
	}
	return nil
}

func checkSkipSamples(sampleID string, skipWordsList []string) bool {
	var skipword bool

	for _, val := range skipWordsList {
		if strings.Contains(strings.ToLower(sampleID), val) {
			skipword = true
			break
		}
	}
	return skipword
}

// Because we stream in the vcf file, we need a way to keep track of
// what columns have the sample ids. We can store the indices in a map
// so that we can get the id later
func mapSamples(header_line string, skipWords []string) map[int]string {
	samplesMap := make(map[int]string)

	split_line := strings.Split(strings.TrimSpace(header_line), "\t")

	for indx, ind_id := range split_line[9:] {
		if checkSkipSamples(ind_id, skipWords) {
			continue
		}
		person_pos := indx + 9
		samplesMap[person_pos] = ind_id
	}

	return samplesMap
}
