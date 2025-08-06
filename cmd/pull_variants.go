package cmd

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/spf13/cobra"
)

var pull_var_cmd = &cobra.Command{
	Use:   "pull-variants",
	Short: "pull variants for the specified region",
	Run:   pull_variants,
}

type VariantInfo struct {
	VariantID   string
	InfoFields  []string
	Calls       string
	Annotations string
}

//	func NewVariantInfo(varID string, variantInfoList []string) *VariantInfo {
//		return &VariantInfo{
//			VariantID:   varID,
//			VariantInfo: variantInfoList,
//			Carrier_map: make(map[string]string),
//		}
//	}
//
// We can parse the genotype calls and determine if there was a no reference call for any of the samples
func parse_genotype_calls(calls []string) bool {
	non_ref_calls := false
	for _, call := range calls {
		switch call {
		case "0/0", "./.", "./0", "0/.", ".":
			continue
		default:
			non_ref_calls = true
		}
	}
	return non_ref_calls
}

func map_header_ids(header_line []string) map[string]int {
	id_mappings := make(map[string]int)

	for indx, id := range header_line {
		id_mappings[id] = indx
	}
	return id_mappings
}

func check_allele_freq(token string, max_freq_threshold float64) (bool, error) {
	maf_field := strings.Split(token, ";")[2]

	maf_values := strings.Split(maf_field, "=")

	for _, maf := range maf_values[1:] {
		// I think the smallest value that a float32 can be is like 1.17e-38 so we should be
		// safe using a 32 bit float because allele frequencies can't get that low in any modern
		// BioBank cohort
		float_val, err := strconv.ParseFloat(maf, 32)
		if err != nil {
			return false, err
		}

		if float_val <= max_freq_threshold {
			return true, nil
		}
	}

	return false, nil
}

func parse_vcf_file(vcf_scanner *bufio.Scanner, maf_cap float64, annotations map[string]string, samples []string, ch chan<- VariantInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	// we also should keep track of lines in the vcf file for error handling
	lineNumber := 0

	// we are going to initialize a hashmap to map the IDs to integers for later on in our loop
	var id_mappings map[string]int
	// now we can parse through the vcf file
	for vcf_scanner.Scan() {

		line := vcf_scanner.Text()

		// increment the line number to the current line we are on
		lineNumber++

		// we can first skip all the unnessecary header lines that have runtime information that we don't need
		if strings.Contains(line, "##") {
			continue
		} else if strings.Contains(line, "#CHROM") {
			// for this line we need to determine the indices that the samples are at in the vcf file line
			split_header := strings.Split(strings.TrimSpace(line), "\t")
			id_mappings = map_header_ids(split_header[9:])
			fmt.Printf("generated id mappings for %d individuals\n", len(id_mappings))
		} else { // If we are not in any of the header lines then we actually want to process the calls
			// We need to make sure the variants are within our region of interest
			split_line := strings.Split(strings.TrimSpace(line), "\t")

			// we also need to get the minor allele freq
			// If there is an error then we can continue in the loop
			if pass_af_threshold, freq_err := check_allele_freq(split_line[7], maf_cap); pass_af_threshold && freq_err == nil {
				// we only need to determine if any of the calls are non variant and then we can return those sites
				if non_ref_call_found := parse_genotype_calls(split_line[9:]); non_ref_call_found {
					// we can build the calls string we need to ensure that the calls are
					// in the same order as the samples with whatever scores we provided
					call_string := strings.Builder{}

					for _, sample_id := range samples {
						// In the id_mapping the indices are start at 0 but in the file the
						// indices for samples will start at 9 so we need to add 9 to the index
						sample_indx := id_mappings[sample_id] + 9
						call_string.WriteString(fmt.Sprintf("%s\t", split_line[sample_indx]))
					}

					// We also need to pull out the annotations for the variant. If the annotation
					// doesn't exist then we can just use a dash string
					anno, ok := annotations[split_line[2]]
					if !ok {
						anno = "-\n"
					}
					variant := VariantInfo{VariantID: split_line[2], InfoFields: split_line[0:9], Calls: call_string.String(), Annotations: anno}
					ch <- variant
				}
			}
		}
		if vcf_scanner.Err() != nil {
			fmt.Printf("Encountered the following error while attempting to read through the vcf file:\n %s\n", vcf_scanner.Err())
		}
	}
	close(ch)
}

func writeToFile(samples []string, annotation_cols []string, writer *bufio.Writer, ch <-chan VariantInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	// we first ned to build the header string. This will have the first 9 fields that are in every
	// vcf file. Then we will add the columns for the sample ids. Then we will add the columns for
	// the annotation fields
	header_str := strings.Builder{}

	header_str.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")

	header_str.WriteString(strings.Join(samples, "\t"))

	header_str.WriteString("\t")

	header_str.WriteString(strings.Join(annotation_cols, "\t"))

	header_str.WriteString("\n")

	_, header_err := writer.WriteString(header_str.String())

	if header_err != nil {
		fmt.Printf("encountered an error while trying to write the header string, %s, to a file. The cause of this could be a bug in the code or unexpected separators in your data. Flushing all of the current data in the writer to the output file but this file is incomplete.\n", header_str.String())
		writer.Flush()
		os.Exit(1)
	}

	// Now we can read through the information in the channel by pulling out 1 variant at a time
	for variant := range ch {
		// now we can build a string for each variant being returned in the analysis
		output_str := strings.Builder{}

		output_str.WriteString(strings.Join(variant.InfoFields, "\t"))
		output_str.WriteString("\t")
		output_str.WriteString(variant.Calls)
		// This code asumes that the variant.Calls ends with a tab separator so
		// therefore we don't need to add any separator between that string and
		// the annotation string. we also assume that the annotation string assumes
		// in a new line
		output_str.WriteString(variant.Annotations)

		_, variant_err := writer.WriteString(output_str.String())

		if variant_err != nil {
			fmt.Printf("encountered an error while trying to write the output variant string, %s, for the variant object, %+v\n. This error could be the result of a bug in the code or an encoding issue within the data. Flushing all current data in the writer but the output file will be incomplete\n", output_str.String(), variant)
			writer.Flush()
			os.Exit(1)
		}
	}
	writer.Flush()
}

func load_column_mappings(line string) map[string]int {
	column_mappings := make(map[string]int)

	column_list := strings.Split(strings.TrimSpace(line), "\t")

	for indx, value := range column_list {
		column_mappings[value] = indx
	}
	return column_mappings
}

func read_annotations(filepath string, cols_to_grab []string) (map[string]string, error) {
	fmt.Printf("Reading in the annotation file: %s\n", filepath)
	annotations := make(map[string]string)

	var err error

	// we need to open the file, pass to a gzip reader and then pass to a bufio Scanner
	anno_fh, anno_err := os.Open(filepath)

	if anno_err != nil {
		err = fmt.Errorf("encountered the following error while attempting to open the annotation file:\n %s", anno_err)
	}

	defer anno_fh.Close()

	gzip_anno_reader, gzip_err := gzip.NewReader(anno_fh)

	if gzip_err != nil {
		err = fmt.Errorf("encountered the following error while trying to decompress the annotation file:\n %s", gzip_err)
	}

	defer gzip_anno_reader.Close()

	buf := make([]byte, 0, 5120*5120)

	anno_scanner := bufio.NewScanner(gzip_anno_reader)

	anno_scanner.Buffer(buf, 5120*5120)

	column_mappings := make(map[string]int)

	for anno_scanner.Scan() {
		cur_line := anno_scanner.Text()
		// we can skip the normal headers in the vcf that are used to
		// keep track of
		if strings.Contains(cur_line, "##") {
			continue
			// We can use the column header line to map the columns that we want to keep
		} else if strings.Contains(cur_line, "#") {
			column_mappings = load_column_mappings(cur_line)
			// otherwise we can process the variants
		} else {
			split_line := strings.Split(strings.TrimSpace(cur_line), "\t")

			// e are going to build the string with all of the columns that we want
			variant_str_builder := strings.Builder{}
			for _, col := range cols_to_grab {
				if value, ok := column_mappings[col]; ok {
					variant_str_builder.WriteString(fmt.Sprintf("\t%s", split_line[value]))
				}
			}

			annotations[split_line[0]] = variant_str_builder.String()
		}
	}
	if anno_scanner.Err() != nil {
		err = fmt.Errorf("encountered the following error while scanner through the annotations file:\n%s", anno_scanner.Err())
	}
	// If there were no annotations loaded into the map then we need to return an error and let the program terminate
	if len(annotations) == 0 {
		err = fmt.Errorf("there were no annotations loading into the internal annotation hashmap after processing the annotations file. This error may could be because the annotation file is empty. but is more likely that the annotation columns that the user desired to keep are not present in the file (Probably due to a spelling error). Please check your annotation file and make sure that the columns you wish to keep are present in the file and spelled the exact same way")
	}
	if len(annotations) < len(cols_to_grab) {
		fmt.Printf("Expected to find %d annotations fields due to the %d keep columns passed to the program. Instead there were only %d fields loaded into the annotations hashmap. There may have been a typo in one of the column names. This is not a fatal issue so the analysis will continue and the output will be missing one of the expected columns. If you wish to end the analysis hit control-C.\n", len(cols_to_grab), len(cols_to_grab), len(annotations))
	}

	fmt.Printf("Read in %d annotations from the file: %s\n", len(annotations), filepath)
	return annotations, err
}

func read_in_samples(samples_filepath string) ([]string, []string) {
	// we are going to return one array of the sample ids and one array of the
	// sample ids with the score appended to the id. This list will be in the
	// same order
	var sample_ids []string
	var formatted_ids []string

	samples_fh, sample_err := os.Open(samples_filepath)

	if sample_err != nil {
		fmt.Printf("Encountered the following error while trying to open the file %s.\n%s\n", samples_filepath, sample_err)
		os.Exit(1)
	}

	defer samples_fh.Close()

	scanner := bufio.NewScanner(samples_fh)

	// this should only be a 2 column file so we should be okay with the standard buffer
	// We are assuming that the first column is the sample id and the second column is the score
	for scanner.Scan() {
		line := scanner.Text()

		split_line := strings.Split(strings.TrimSpace(line), "\t")

		if len(split_line) == 1 {
			sample_ids = append(sample_ids, split_line[0])
			formatted_ids = append(formatted_ids, split_line[0])
		} else {
			sample_ids = append(sample_ids, split_line[0])

			if dot_indx := strings.Index(split_line[1], "."); dot_indx != -1 {
				trimmed_score := split_line[1][0 : dot_indx+3]
				formatted_id := fmt.Sprintf("%s_%s", split_line[0], trimmed_score)
				formatted_ids = append(formatted_ids, formatted_id)
			} else {
				formatted_id := fmt.Sprintf("%s_%s", split_line[0], split_line[1])
				formatted_ids = append(formatted_ids, formatted_id)
			}
		}
	}
	if scanner.Err() != nil {
		fmt.Printf("Encountered an error while scanning through the samples file:\n %s\n.", scanner.Err())
	}

	fmt.Printf("Read in %d samples and formatted %d IDs\n", len(sample_ids), len(formatted_ids))

	return sample_ids, formatted_ids
}

func pull_variants(cmd *cobra.Command, args []string) {
	start_time := time.Now()

	fmt.Printf("began the analysis at: %s\n", start_time.String())
	// parse all the arguments needs for this command
	anno_file, _ := cmd.Flags().GetString("anno-file")
	cols_to_keep, _ := cmd.Flags().GetString("keep-cols")
	samples_filepath, _ := cmd.Flags().GetString("samples")
	output_file, _ := cmd.Flags().GetString("output")
	maf_cap, _ := cmd.Flags().GetFloat64("maf-threshold")
	// log_filepath, _ := cmd.Flags().GetString("log-filepath")

	// read in the annotations into a dictionary

	anno_cols_to_keep := strings.Split(cols_to_keep, ",")

	anno_map, anno_err := read_annotations(anno_file, anno_cols_to_keep)
	if anno_err != nil {
		fmt.Printf("Encountered the following error while trying to read in the annotations.\n %s\n", anno_err)
		os.Exit(1)
	}

	// we also need to read in the samples file. We are going to return 2 values. One will
	// be the list of ids as we encounter them in the file. The other will be the list of
	// ids with the phers score appended
	samples, formatted_samples := read_in_samples(samples_filepath)

	// lets read from stdin. We need to increase the buffer because the default buffer is too small for our files
	buf := make([]byte, 0, 5120*5120)

	buffered_vcf := bufio.NewScanner(os.Stdin)

	buffered_vcf.Buffer(buf, 5120*5120)
	// We also need to open the output file for writing
	output_fh, output_err := os.Create(output_file)

	if output_err != nil {
		fmt.Printf("There was an issue trying to create the output file: %s\n", output_file)
		os.Exit(1)
	}

	defer output_fh.Close()

	writer := bufio.NewWriter(output_fh)

	// lets create a channel and a waitgroup so we can have the parsing vcf in one goroutine and the writing in another goroutine
	ch := make(chan VariantInfo)
	var wg sync.WaitGroup

	wg.Add(1)
	// now we can parse the vcf file
	go parse_vcf_file(buffered_vcf, maf_cap, anno_map, samples, ch, &wg)

	wg.Add(1)

	go writeToFile(formatted_samples, anno_cols_to_keep, writer, ch, &wg)

	wg.Wait()

	end_time := time.Now()

	fmt.Printf("finished analysis at: %s\n", end_time.String())

	duration := end_time.Sub(start_time)

	fmt.Printf("total analysis time: %s\n", duration.String())
}

// I am assuming that the user is using bcftools to stream data into this program. Therefore
// we only need to read from the stdin stream and don't nedd them to provide the vcf file as
// input
func init() {
	RootCmd.AddCommand(pull_var_cmd)
	pull_var_cmd.Flags().StringP("anno-file", "a", "", "Filepath to write the output file to.")
	pull_var_cmd.Flags().StringP("keep-cols", "c", "", "Columns in the annotation file to keep.")
	pull_var_cmd.Flags().StringP("samples", "s", "", "file where the first column contains the samples for ")
	pull_var_cmd.Flags().StringP("output", "o", "test_output.txt", "Filepath to write the output file to.")
	pull_var_cmd.Flags().String("log-filepath", "text.log", "Filepath to write the log file to.")
	pull_var_cmd.Flags().Float64("maf-threshold", 0.1, "Minor allele frequency cap to filter output to only variants less than this value")
}
