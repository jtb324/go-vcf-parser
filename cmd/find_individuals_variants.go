package cmd

import (
	"bufio"
	"fmt"
	"os"
	"slices"
	"strings"
	"time"

	"github.com/spf13/cobra"
)

var sample_var_cmds = &cobra.Command{
	Use:   "view-sample-variants",
	Short: "grab the variants that samples of interest have. This command uses the output from the pull-variants command",
	Run:   find_sample_variants,
}

type SampleInfo struct {
	Score    string
	Variants []string
}

type SampleID struct {
	Index    int
	SampleID string
	Score    string
}

func read_samples_file(samples_filepath string) ([]string, []error) {
	fmt.Printf("Reading in all of the desired samples from the file %s\n", samples_filepath)
	var errors []error
	var samples []string

	samples_fh, samples_err := os.Open(samples_filepath)
	if samples_err != nil {
		errors = append(errors, fmt.Errorf("failed to open the file, %s. The following error was encountered, %s", samples_filepath, samples_err))
	} else {
		sample_scanner := bufio.NewScanner(samples_fh)
		for sample_scanner.Scan() {
			line := sample_scanner.Text()
			if strings.Contains(strings.ToLower(line), "grid") || strings.Contains(strings.ToLower(line), "IID") {
				continue
			}
			split_line := strings.Split(strings.TrimSpace(line), "\t")

			samples = append(samples, split_line[0])
		}
		if sample_scanner.Err() != nil {
			errors = append(errors, fmt.Errorf("encountered the following error while scanning through the samples file. %s", sample_scanner.Err()))
		}
	}

	// If no samples were read in the we should give the user an error
	if len(samples) == 0 {
		errors = append(errors, fmt.Errorf("unable to read any of the samples from the file, %s. Please ensure that the file is tab separated and that each row is 1 individual", samples_filepath))
	} else {
		fmt.Printf("Read %d samples in from the file, %s", len(samples), samples_filepath)
	}
	return samples, errors
}

func parse_header(header_line string, samples []string) []SampleID {
	var header_map []SampleID

	split_header := strings.Split(strings.TrimSpace(header_line), "\t")

	for indx, val := range split_header {
		// Sometimes the id will have a score (either PheRS or case/control status) appended to the end. We can split the string to get this value
		split_id := strings.Split(val, "_")

		// we can check if the sample id (split_id[0]) is in our samples array. If not then we skip that position
		if !slices.Contains(samples, split_id[0]) {
			continue
		}
		if len(split_id) == 2 {
			header_map = append(header_map, SampleID{Index: indx, SampleID: split_id[0], Score: split_id[1]})
		} else {
			header_map = append(header_map, SampleID{Index: indx, SampleID: split_id[0], Score: ""})
		}
	}
	return header_map
}

func check_for_alt_call(call string, reference_call_set map[string]bool) bool {
	_, ok := reference_call_set[call]
	// If the call is an alt call then the dictionary will return
	// false because it only contains keys for reference calls.
	// We can negate the false value so that it equals true
	// indicating that we did find an alternate allele
	return !ok
}

func parse_calls(samples_scanner *bufio.Scanner, samples []string) (map[string]SampleInfo, []error) {
	sampleInfo := make(map[string]SampleInfo) // This will be our return value
	var errors []error
	var header_map []SampleID
	// We also need to generate the set of reference calls so that we can compare our calls for that
	reference_call_strs := generate_reference_set()
	// This file has a header line so we first need to read in the indices for each column
	for samples_scanner.Scan() {
		line := samples_scanner.Text()
		// We assume the header line contains the phrase #CHROM because this is the output of the other program
		if strings.Contains(line, "#CHROM") {
			header_map = parse_header(line, samples)
		} else {
			split_line := strings.Split(strings.TrimSpace(line), "\t")
			for _, individual := range header_map {
				call := split_line[individual.Index]
				if check_for_alt_call(call, reference_call_strs) {
					variantStr := fmt.Sprintf("%s:%s", split_line[1], call)
					// We need to keep track of the variants that an individual has. The variant is the first value in the split_line array.
					if sampleStruct, ok := sampleInfo[individual.SampleID]; ok {
						sampleStruct.Variants = append(sampleStruct.Variants, variantStr)
					} else {
						variantList := []string{variantStr}
						sampleInfo[individual.SampleID] = SampleInfo{Score: individual.Score, Variants: variantList}
					}
				}
			}
		}
	}
	if samples_scanner.Err() != nil {
		errors = append(errors, fmt.Errorf("encountered the following error while trying to scan through the calls file:  %s", samples_scanner.Err()))
	}

	return sampleInfo, errors
}

func write_variants(writer *bufio.Writer, sample_variants map[string]SampleInfo) {
	// lets build the header line

	header_str := "SAMPLE\tSCORE\tVARIANTS\n"

	writer.WriteString(header_str)

	for sample_id, sampleInfoObj := range sample_variants {
		sample_str := strings.Builder{}

		sample_str.WriteString(sample_id)
		// We can build the rest of the string appending the Score if there is one and the variants
		if sampleInfoObj.Score == "" {
			sample_str.WriteString(fmt.Sprintf("\t-\t%s", strings.Join(sampleInfoObj.Variants, ",")))
		} else {
			sample_str.WriteString(fmt.Sprintf("\t%s\t%s", sampleInfoObj.Score, strings.Join(sampleInfoObj.Variants, ",")))
		}
		writer.WriteString(sample_str.String())
	}
}

func find_sample_variants(cmd *cobra.Command, args []string) {
	start_time := time.Now()

	fmt.Printf("began the analysis at: %s\n", start_time.Format("2006-01-02@15:04:05"))
	// read in the appropriate CLI flags
	calls_file, _ := cmd.Flags().GetString("calls-file")
	samples_list, _ := cmd.Flags().GetString("samples-list")
	samples_filepath, _ := cmd.Flags().GetString("samples-file")
	output_filepath, _ := cmd.Flags().GetString("output")

	var samples []string
	var sample_file_err []error
	if samples_list != "" {
		samples = strings.Split(samples_list, ",")
	} else if samples_filepath != "" {
		// process the samples file
		samples, sample_file_err = read_samples_file(samples_filepath)
		if sample_file_err != nil {
			fmt.Printf("Encountered the following errors while trying to read in samples from the file %s\n", samples_filepath)
			for msg_indx, msg := range sample_file_err {
				fmt.Printf("Error Msg %d:\n %s", msg_indx, msg)
			}
			os.Exit(1)
		}
	}
	// now we can parse through the output file for variants of interest
	calls_fh, calls_err := os.Open(calls_file)

	if calls_err != nil {
		fmt.Printf("Encountered the following error while trying to open the file, %s\n%s\n", calls_file, calls_err)
		os.Exit(1)
	}
	// Create the scanner to read the calls file with a custom buffer
	buf := make([]byte, 0, 1024*1024)

	calls_scanner := bufio.NewScanner(calls_fh)

	calls_scanner.Buffer(buf, 1024*1024)

	sample_variants, errs := parse_calls(calls_scanner, samples)

	fmt.Printf("Identified variants for %d samples\n", len(sample_variants))

	if errs != nil {
		fmt.Println("Encountered the following errors while scanning through the calls file to identify the variants that each sample has")
		for indx, err_msg := range errs {
			fmt.Printf("Error Msg %d:\n %s\n", indx, err_msg)
		}
		os.Exit(1)
	}

	output_fh, output_err := os.Open(output_filepath)

	if output_err != nil {
		fmt.Printf("Encountered the following error while trying to open the output file, %s.\n %s\n", output_filepath, output_err)
		os.Exit(1)
	}

	writer := bufio.NewWriter(output_fh)
	fmt.Printf("Writing output to the file: %s\n", output_filepath)
	write_variants(writer, sample_variants)
}

// I am assuming that the user is using bcftools to stream data into this program. Therefore
// we only need to read from the stdin stream and don't nedd them to provide the vcf file as
// input
func init() {
	RootCmd.AddCommand(sample_var_cmds)
	sample_var_cmds.Flags().StringP("calls-file", "c", "", "output tab separated file from the pull-variants command. Each row should be a variant and each column will be an individuals with the phers score (if the score was provided in that command)")
	sample_var_cmds.Flags().StringP("samples-list", "s", "", "list of sample ids to find all the variants for. This list should be comma separated with now spaces inbetween ids.")
	sample_var_cmds.Flags().StringP("samples-file", "S", "", "filepath to a tab separated text file that has the samples we wish to keep. The first column is expected to be a list of grids and the file may or may not have a header")
	sample_var_cmds.Flags().StringP("output", "o", "test_output.txt", "Filepath to write the output file to.")
	sample_var_cmds.Flags().String("log-filepath", "text.log", "Filepath to write the log file to.")
}
