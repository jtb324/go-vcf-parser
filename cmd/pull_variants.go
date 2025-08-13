package cmd

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"time"

	gzip "github.com/klauspost/pgzip"
	"github.com/spf13/cobra"
)

var pull_var_cmd = &cobra.Command{
	Use:   "pull-variants",
	Short: "pull variants for the specified region",
	Run:   pull_variants,
}

type VariantAnnotations map[string]*strings.Builder

type VariantInfo struct {
	VariantID   string
	InfoFields  []string
	Calls       string
	Annotations VariantAnnotations
}

func generate_reference_set() map[string]bool {
	ref_call := map[string]bool{
		"0/0": true,
		"./.": true,
		"0/.": true,
		"./0": true,
		".":   true,
	}

	return ref_call
}

// We can parse the genotype calls and determine if there was a no reference call for any of the samples
// The ref calls is a dictionary where the keys are the different reference genotype calls. If we try to use the sample calls as errors and we find one that is not a valid key then it indicates that we have found a non-reference call
func parse_genotype_calls(calls []string, ref_calls map[string]bool) bool {
	non_ref_calls := false
	for _, call := range calls {
		if _, ok := ref_calls[call]; !ok {
			non_ref_calls = true
			break
		}
	}
	return non_ref_calls
}

func map_header_ids(samples []string) map[string]int {
	id_mappings := make(map[string]int)

	for indx, id := range samples {
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

func process_header_ids(vcf_scanner *bufio.Scanner, pheno_map map[string]string) ([]string, string, error) {
	// We need to return a list of the samples. This value will be used while parsing the vcf file sequencing calls.
	var samples []string
	// create the sample string builder so that we can add ids as we process them. This string will be used when writting the output
	sample_str := strings.Builder{}

	var err error
	samples_count := 0 // We also are going to keep counts of the number of samples so that we can report that back to the user

	// starting a counter for the line number which can be used in error messages
	line_number := 0
Scanner: // we can create a label for the outer scanner loop so that we can selectively break all the way out
	for vcf_scanner.Scan() {
		line := vcf_scanner.Text()

		line_number++

		if strings.Contains(line, "##") {
			continue
		} else if strings.Contains(line, "#CHROM") {
			split_header := strings.Split(strings.TrimSpace(line), "\t")
			// we can now set the samples
			samples = split_header[9:]
			for _, id := range split_header[9:] { // sample IDs start at the 9 index in the vcf file. This is standard format
				if value, ok := pheno_map[id]; ok {
					sample_str.WriteString(fmt.Sprintf("%s_%s\t", id, value))
					samples_count++
				} else {
					err = fmt.Errorf("the id %s had no phenotype information meaning that it was not present in the phenotype file but it is present in the header of the VCF file that is being streamed in. This error may be the result of providing an incorrect version of either the phenotype file to the program or the samples file used to filter from bcftools. Please rectify this two files so that the samples file either has the same individuals as the phenotype file or it is a subset of the individuals in the phenotype file. Program will now terminate", id)
					break Scanner // Break out of the whole scanner loop
				}
			}
			fmt.Printf("processed the header line for the provided vcf file and identified %d samples in the header\n", samples_count)
			break Scanner // we don't want to process more of the scanner so we can break after processing the header
		} else { // If we are not in any of the header lines then we actually want to process the calls
			err = fmt.Errorf("there was no header detected in the initial lines of the vcf file. This program expects for there to be a header line that begins with a single #CHROM in order to identify the sample ids. This lack of header may be the result of streaming the data using bcftools with the -H flag to remove the header. Please remove this flag, and rerun the program")
			break Scanner // If there is no header then we don't need to process the file. We just want to return an error and then terminate the program
		}
	}
	if vcf_scanner.Err() != nil {
		err = fmt.Errorf("encountered the following error on line %d while trying to scan through the header of the vcf file for sample ids: %s", line_number, vcf_scanner.Err())
	}
	// The final sample_str will end in a tab separator. This needs to be kept in mind when writing the string to a file
	return samples, sample_str.String(), err
}

func parse_vcf_file(vcf_scanner *bufio.Scanner, maf_cap float64, annotations map[string]VariantAnnotations, samples []string, sample_indices map[string]int, ch chan<- VariantInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	// Lets create the reference genotype map
	reference_calls := generate_reference_set()
	// now we can parse through the vcf file. We don't have to account for the header lines
	// because we have a separator function handling this before the go routines
	for vcf_scanner.Scan() {

		line := vcf_scanner.Text()

		// we can first skip all the unnessecary header lines that have runtime information that we don't need
		// We need to make sure the variants are within our region of interest
		split_line := strings.Split(strings.TrimSpace(line), "\t")

		// we also need to get the minor allele freq
		// If there is an error then we can continue in the loop
		if pass_af_threshold, freq_err := check_allele_freq(split_line[7], maf_cap); pass_af_threshold && freq_err == nil {
			// we only need to determine if any of the calls are non variant and then we can return those sites
			if non_ref_call_found := parse_genotype_calls(split_line[9:], reference_calls); non_ref_call_found {
				// we can build the calls string we need to ensure that the calls are
				// in the same order as the samples with whatever scores we provided
				call_string := strings.Builder{}

				for _, sample_id := range samples {
					// In the id_mapping the indices are start at 0 but in the file the
					// indices for samples will start at 9 so we need to add 9 to the index
					sample_indx := sample_indices[sample_id] + 9
					call_string.WriteString(fmt.Sprintf("\t%s", split_line[sample_indx]))
				}

				// We also need to pull out the annotations for the variant. If the annotation
				// doesn't exist then we can just use an empty string. The ok returns true if
				// the value is in the dictionary and false if it is not.
				anno, ok := annotations[split_line[2]]
				if !ok {
					anno = nil
				}
				variant := VariantInfo{VariantID: split_line[2], InfoFields: split_line[0:9], Calls: call_string.String(), Annotations: anno}
				ch <- variant
			}
		}
		if vcf_scanner.Err() != nil {
			fmt.Printf("Encountered the following error while attempting to read through the vcf file:\n %s\n", vcf_scanner.Err())
		}
	}
	close(ch)
}

// parse the VariantAnnotations
func generate_annotation_str(variant_annos VariantAnnotations, anno_cols []string) string {
	annotation_str := strings.Builder{}
	for _, col := range anno_cols {
		if value, ok := variant_annos[col]; ok {
			formatted_val := fmt.Sprintf("\t%s", value.String())
			annotation_str.WriteString(formatted_val)
		}
	}
	return annotation_str.String()
}

func writeToFile(samples string, annotation_cols []string, writer *bufio.Writer, ch <-chan VariantInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	// counter to record how many variants were written to a file
	variants_written := 0
	// we first ned to build the header string. This will have the first 9 fields that are in every
	// vcf file. Then we will add the columns for the sample ids. Then we will add the columns for
	// the annotation fields
	header_str := strings.Builder{}

	header_str.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")

	header_str.WriteString(samples)

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
		// WE first join initial 9 fields from the vcf file that we stored in the variant.InfoFields attribute
		output_str.WriteString(strings.Join(variant.InfoFields, "\t"))
		// next we can append the calls to this string. This calls string starts with a tab character
		output_str.WriteString(variant.Calls)
		// This code asumes that the variant.Calls ends with a tab separator so
		// therefore we don't need to add any separator between that string and
		// If the annotation string is empty then there were no annotations for the specific variant
		// and we have to create the annotation string by just creating '-' for each column
		if variant.Annotations == nil {
			for range annotation_cols {
				output_str.WriteString("\t-")
			}
			output_str.WriteString("\n")
		} else {
			anno_str := generate_annotation_str(variant.Annotations, annotation_cols)
			output_str.WriteString(fmt.Sprintf("%s\n", anno_str))
		}

		_, variant_err := writer.WriteString(output_str.String())

		if variant_err != nil {
			fmt.Printf("encountered an error while trying to write the output variant string, %s, for the variant object, %+v\n. This error could be the result of a bug in the code or an encoding issue within the data. Flushing all current data in the writer but the output file will be incomplete\n", output_str.String(), variant)
			writer.Flush()
			os.Exit(1)
		}
		// increment the variants_written counter to represent that we have written another variant to file
		variants_written++
	}
	writer.Flush()
	fmt.Printf("Recorded information for %d variant(s)\n", variants_written)
}

func load_column_mappings(line string) map[string]int {
	column_mappings := make(map[string]int)

	column_list := strings.Split(strings.TrimSpace(line), "\t")

	for indx, value := range column_list {
		column_mappings[value] = indx
	}
	return column_mappings
}

func check_region(anno_pos string, start int, end int) (bool, []error) {
	re := regexp.MustCompile(`[:|-]`)
	split_pos := re.Split(anno_pos, -1)

	// This makes the assumption that the position either has the from chr:pos or pos or chr:pos1-pos2.
	// If the position is only "pos" then the split will produce an array of 1 value. If there is only
	// a : then it will have 2 values. If it has both a : and - then the resulting array will have 3
	// values. If the length of the resulting array is > 2 then we just need to pull the second value. If the length is 3 then we need to set the start and end
	var start_pos_str string
	var end_pos_str string
	var conversion_err []error

	if len(split_pos) == 1 {
		start_pos_str = split_pos[0]
	} else if len(split_pos) == 2 {
		start_pos_str = split_pos[1]
	} else {
		start_pos_str = split_pos[1]
		end_pos_str = split_pos[2]
	}

	start_pos, first_conv_err := strconv.Atoi(start_pos_str)

	if first_conv_err != nil {
		conversion_err = append(conversion_err, fmt.Errorf("encountered the following error while converting the starting position of the string, %s\n. %s", anno_pos, first_conv_err))
	}

	// end_pos_str may be an empty string at this point which will cause the conversion to fail. The
	// empty string occurs because there is not always a second position in the anno_pos string. If
	// failure occurs because of this we don't actually care and we can still use the empty value of
	// end_pos (0) in our logical statement to make the second test fail and shortcircuit
	end_pos, second_conv_err := strconv.Atoi(end_pos_str)

	if end_pos_str != "" && second_conv_err != nil {
		conversion_err = append(conversion_err, fmt.Errorf("enocuntered the following error while converting the ending position of the string %s\n. %s", anno_pos, second_conv_err))
	}
	// Here we are going to check if the starting region falls within our desired region or if the end end position falls within our starting region
	return start <= start_pos && start_pos <= end || (end_pos != 0 && start <= end_pos && end_pos <= end), conversion_err
}

func read_annotations(filepath string, cols_to_grab []string, region Region) (map[string]VariantAnnotations, error) {
	fmt.Printf("Reading in the annotation file: %s\n", filepath)
	fmt.Printf("Collecting annotations only for sites overlapping this region: %s:%d-%d\n", region.chrom, region.start, region.end)
	annotations := make(map[string]VariantAnnotations)

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

	buf := make([]byte, 0, 7168*7168)

	// anno_scanner := bufio.NewScanner(gzip_anno_reader)
	anno_scanner := bufio.NewScanner(gzip_anno_reader)

	anno_scanner.Buffer(buf, 7168*7168)

	column_mappings := make(map[string]int)
Main_Loop:
	for anno_scanner.Scan() {
		cur_line := anno_scanner.Text()
		// we can skip the normal headers in the vcf that are used to
		// log information from prior transformations done on the data
		if strings.Contains(cur_line, "##") {
			continue
			// We can use the column header line to map the columns that we want to keep
		} else if strings.Contains(cur_line, "#") {
			column_mappings = load_column_mappings(cur_line)
			// otherwise we can process the variants
		} else {
			// Once we are past all of the header lines then we can pull information for each variant.
			// Sometimes variants also have multiple transcripts and therefore show up on multiple rows.
			// We have to handle this by aggregating together the different information
			// we can use a string builder to keep track of the annotation and separate the different values by a comma
			split_line := strings.Split(strings.TrimSpace(cur_line), "\t")
			// first lets see if this annotation is even in the right position. If it is not in the right position then we can just continue the loop
			if in_region, ok := check_region(split_line[1], region.start, region.end); !in_region && ok == nil {
				continue Main_Loop
			} else if ok != nil {
				fmt.Printf("Encountered an issue while checking if the variant %s was in the search region of %d-%d\n %s\n Skipping this variant and proceeding to the next one", split_line[1], region.start, region.end, ok)
			}
			// we can check if there is already an annotation created for the variant and add things to it. Otherwise we can just
			variant_annotations := annotations[split_line[0]]
			// if the anotation is present then we can iterate over the columns and update the string.builder for each appropriate columns
			if variant_annotations != nil {
				for _, col := range cols_to_grab {
					if value, ok := column_mappings[col]; ok {
						value_str := fmt.Sprintf(";%s", split_line[value])
						variant_annotations[col].WriteString(value_str)
					}
				}
				// otherwise we have to create a new map that will have a key for each column in the
				// analysis. We can then iterate over each column and append information to the string.Builder for that key
			} else {
				variant_annos := make(VariantAnnotations)
				for _, col := range cols_to_grab {
					col_values := strings.Builder{}
					if value, ok := column_mappings[col]; ok {
						col_values.WriteString(split_line[value])
						variant_annos[col] = &col_values
					}
				}
				annotations[split_line[0]] = variant_annos
			}
		}
	}
	if anno_scanner.Err() != nil {
		err = fmt.Errorf("encountered the following error while scanner through the annotations file:\n%s", anno_scanner.Err())
	}
	// If there were no annotations loaded into the map then we need to return an error and let the program terminate
	if len(annotations) == 0 {
		err = fmt.Errorf("there were no annotations loading into the internal annotation hashmap after processing the annotations file. This error may could be because the annotation file is empty. but is more likely that the annotation columns that the user desired to keep are not present in the file (Probably due to a spelling error). Please check your annotation file and make sure that the columns you wish to keep are present in the file and spelled the exact same way")
	}
	// if len(annotations) < len(cols_to_grab) { // TODO: This error message is actually not right. Currently it gets set off it there are no annotations found.
	//
	// 	fmt.Printf("Expected to find %d annotations fields due to the %d keep columns passed to the program. Instead there were only %d fields loaded into the annotations hashmap. There may have been a typo in one of the column names. This is not a fatal issue so the analysis will continue and the output will be missing one of the expected columns. If you wish to end the analysis hit control-C.\n", len(cols_to_grab), len(cols_to_grab), len(annotations))
	// }

	fmt.Printf("Read in %d annotations from the file: %s\n", len(annotations), filepath)
	return annotations, err
}

func read_in_samples(samples_filepath string) map[string]string {
	// we are going to return one array of the sample ids and one array of the
	// sample ids with the score appended to the id. This list will be in the
	// same order
	sample_ids := make(map[string]string)

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
			sample_ids[split_line[0]] = ""
		} else {
			if dot_indx := strings.Index(split_line[1], "."); dot_indx != -1 {
				trimmed_score := split_line[1][0 : dot_indx+3]
				sample_ids[split_line[0]] = trimmed_score
			} else {
				sample_ids[split_line[0]] = split_line[1]
			}
		}
	}
	if scanner.Err() != nil {
		fmt.Printf("Encountered an error while scanning through the samples file:\n %s\n.", scanner.Err())
	}

	fmt.Printf("Read in %d samples from the file: %s\n", len(sample_ids), samples_filepath)

	return sample_ids
}

type Region struct {
	chrom string
	start int
	end   int
}

func parse_region(region_str string) (Region, []error) {
	re := regexp.MustCompile(`[:|-]`)
	region_split := re.Split(region_str, -1)

	var err []error
	var region Region

	if len(region_split) == 1 {
		err = append(err, fmt.Errorf("failed to split the region string. Make sure that the region string is of the form chrX:start-end"))
	} else {

		start_int, start_err := strconv.Atoi(region_split[1])

		if start_err != nil {
			err = append(err, fmt.Errorf("encountered the following error when we tried to convert the starting position of the region string %s to an integer: %s", region_str, start_err))
		}

		end_int, end_err := strconv.Atoi(region_split[2])

		if end_err != nil {
			err = append(err, fmt.Errorf("encountered the following error when we tried to convert the ending position of the region string %s to an integer: %s", region_str, end_err))
		}
		// We do need to make sure that the end point is not smaller than the start point because that will mess many things up
		if start_int >= end_int {
			err = append(err, fmt.Errorf("the parsed end point is smaller than the starting point of the region. This suitation will result in no annotations being loaded from the annotation file later on. This issue may mean that there is a typo in the region flag. Please check this flag and make sure that the end position is greater than the start position"))
		}
		region = Region{chrom: region_split[0], start: start_int, end: end_int}
	}
	return region, err
}

func pull_variants(cmd *cobra.Command, args []string) {
	start_time := time.Now()

	fmt.Printf("began the analysis at: %s\n", start_time.Format("2006-01-02@15:04:05"))
	// parse all the arguments needs for this command
	anno_file, _ := cmd.Flags().GetString("anno-file")
	cols_to_keep, _ := cmd.Flags().GetString("keep-cols")
	samples_filepath, _ := cmd.Flags().GetString("samples")
	output_file, _ := cmd.Flags().GetString("output")
	maf_cap, _ := cmd.Flags().GetFloat64("maf-threshold")
	region, _ := cmd.Flags().GetString("region")
	// log_filepath, _ := cmd.Flags().GetString("log-filepath")
	// lets parse the region
	parsed_region, region_err := parse_region(region)

	if region_err != nil {
		fmt.Printf("Encountered the following errors while trying to parse the region value: \n")
		for _, msg := range region_err {
			fmt.Println(msg)
		}
		// These issues are all worth terminating the program
		os.Exit(1)
	}
	// read in the annotations into a dictionary

	anno_cols_to_keep := strings.Split(cols_to_keep, ",")

	anno_map, anno_err := read_annotations(anno_file, anno_cols_to_keep, parsed_region)
	if anno_err != nil {
		fmt.Printf("Encountered the following error while trying to read in the annotations.\n %s\n", anno_err)
		os.Exit(1)
	}

	// we also need to read in the samples file. We are going to return 2 values. One will
	// be the list of ids as we encounter them in the file. The other will be the list of
	// ids with the phers score appended
	sample_phenos := read_in_samples(samples_filepath)

	// lets read from stdin. We need to increase the buffer because the default buffer is too small for our files
	buf := make([]byte, 0, 5120*5120)

	buffered_vcf := bufio.NewScanner(os.Stdin)

	buffered_vcf.Buffer(buf, 5120*5120)

	// We need to process the header row first. Ids in the sample string are in the same
	// order as the samples but they have the phenotype information added to the string
	// formatted as "_score"
	samples, sample_str, header_err := process_header_ids(buffered_vcf, sample_phenos)

	if header_err != nil {
		fmt.Printf("%s\nTerminating programming...\n", header_err)
		os.Exit(1)
	}
	// we then nedd to use the samples list and map this values to an index because
	// this is the order they will be in the vcf stream
	samples_indices := map_header_ids(samples)
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
	go parse_vcf_file(buffered_vcf, maf_cap, anno_map, samples, samples_indices, ch, &wg)

	wg.Add(1)

	go writeToFile(sample_str, anno_cols_to_keep, writer, ch, &wg)

	wg.Wait()

	end_time := time.Now()

	fmt.Printf("finished analysis at: %s\n", end_time.Format("2006-01-02@15:04:05"))

	duration := end_time.Sub(start_time)

	fmt.Printf("total analysis time: %s\n", duration.String())
}

// I am assuming that the user is using bcftools to stream data into this program. Therefore
// we only need to read from the stdin stream and don't nedd them to provide the vcf file as
// input
func init() {
	RootCmd.AddCommand(pull_var_cmd)
	pull_var_cmd.Flags().StringP("anno-file", "a", "", "Filepath to an annotation file (currently only supports VEP so that there is a canocial column that we can use to avoid duplicates and only look at the cannocial transcript).")
	pull_var_cmd.Flags().StringP("keep-cols", "c", "", "Columns in the annotation file to keep.")
	pull_var_cmd.Flags().StringP("samples", "s", "", "file where the first column contains the samples for ")
	pull_var_cmd.Flags().StringP("output", "o", "test_output.txt", "Filepath to write the output file to.")
	pull_var_cmd.Flags().StringP("region", "r", "", "Region of the chromosome that we are interested in searching for. This region should have the form chrX:start-end. We will use this region to filter which annotations we save in memory")
	pull_var_cmd.Flags().String("log-filepath", "text.log", "Filepath to write the log file to.")
	pull_var_cmd.Flags().Float64("maf-threshold", 0.1, "Minor allele frequency cap to filter output to only variants less than this value")
}
