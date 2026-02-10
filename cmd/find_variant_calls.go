package cmd

import (
	"bufio"
	"fmt"
	"go-phers-parser/internal/files"
	"maps"
	"os"
	"slices"
	"strings"
)

func check_alt_call(call string, reference_call_set map[string]bool) bool {
	_, call_is_ref := reference_call_set[call] // line checks to see if our value is one of the reference calls

	return !call_is_ref
}

type Result struct {
	Variants []VariantCalls
	Errors   []error
	Samples  map[string]bool
}

func (result *Result) generate_sample_list() []string {
	return slices.Collect(maps.Keys(result.Samples))
}

type VariantCalls struct {
	VariantInfo     []string
	VariantCarriers map[string]string
	GenotypeCounts  map[string]int
}

func update_genotype_count(call string, genotype_counts map[string]int) {
	switch call {
	case "0/0":
		genotype_counts["homo_ref"]++
	case "0/1", "1/0":
		genotype_counts["het"]++
	case "1/1":
		genotype_counts["homo_alt"]++
	case "./.", ".":
		genotype_counts["no_calls"]++
	default:
		genotype_counts["other"]++
	}
}

func process_variant_stream(streamReader *files.VCFReader, resultsObj *Result) error {
	for streamReader.FileScanner.Scan() {

		// We can initialize the variantCalls object with a dictionary for the genotype counts.
		// This structure will help us while writing later
		variantCallsObj := VariantCalls{
			VariantCarriers: make(map[string]string),
			GenotypeCounts: map[string]int{
				"homo_alt": 0,
				"homo_ref": 0,
				"het":      0,
				"no_calls": 0,
				"other":    0,
			},
		}

		line := streamReader.FileScanner.Text()
		split_line := strings.Split(strings.TrimSpace(line), "\t")

		// We can add the variant string here
		variantCallsObj.VariantInfo = split_line[0:3]

		// We will need to generate the reference calls for comparison
		ref_call_set := generate_reference_set()
		// We can iterate over each call
		for indx, calls := range split_line[9:] {
			indx = indx + 9
			// There may be some indices that are missing if there are samples we want to skip.
			// We will need to check and make sure the key exist and only proceed if it does
			if id, ok := streamReader.SampleMapping[indx]; ok {
				if check_alt_call(calls, ref_call_set) {
					// We can add the id and the call to the carriers map
					variantCallsObj.VariantCarriers[id] = calls
					// Then we can also save the carrier ids we found. We will use
					// this list to create the header for the output file later
					resultsObj.Samples[id] = true // This is how you use a set in Go. Its the same as a map
				}
				update_genotype_count(calls, variantCallsObj.GenotypeCounts)
			}
		}
		fmt.Printf("Identified %d individuals who were either heterozygous or homozygous alt for the variant %s\n", len(variantCallsObj.VariantCarriers), variantCallsObj.VariantInfo[2])
		resultsObj.Variants = append(resultsObj.Variants, variantCallsObj)
	}
	if streamReader.FileScanner.Err() != nil {
		return streamReader.FileScanner.Err()
	}
	return nil
}

func writer(writer *bufio.Writer, results Result) {
	// get a list of all the samples we need to put in the header
	sample_list := results.generate_sample_list()
	// Create the header string
	header_str := strings.Builder{}
	header_str.WriteString("CHROM\tPOS\tID\tHOMO_REF_COUNT\tHET_COUNT\tHOMO_ALT_COUNT\tNO_CALL_COUNT\tOTHER_CALL_COUNT\t")
	header_str.WriteString(fmt.Sprintf("%s\n", strings.Join(sample_list, "\t")))

	writer.WriteString(header_str.String())
	// Now create the output string
	for _, variant := range results.Variants {
		row_str := strings.Builder{}
		row_str.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%d", strings.Join(variant.VariantInfo, "\t"), variant.GenotypeCounts["homo_ref"], variant.GenotypeCounts["het"], variant.GenotypeCounts["homo_alt"], variant.GenotypeCounts["no_calls"], variant.GenotypeCounts["other"]))
		for sampleID := range results.Samples {
			sample_call, ok := variant.VariantCarriers[sampleID]

			var output_str string
			if ok {
				output_str = fmt.Sprintf("\t%s:%s", sampleID, sample_call)
			} else {
				output_str = "\t-"
			}
			row_str.WriteString(output_str)
		}
		row_str.WriteString("\n")
		writer.WriteString(row_str.String())
	}
	writer.Flush()
}

// This function is used to find all the individuals with variant calls for a site of interest.
// It expects to have input streamed in from bcftools
func FindAllCarrierCalls(output_filepath string, buffersize int, exclusion_substring string) {

	// we need to create the reader
	vcfStreamer := files.MakeStreamReader(buffersize)

	// We need to add the sample-exclusion-string
	vcfStreamer.SampleExclusions = strings.Split(exclusion_substring, ",")

	// We need to early terminate if there was an error while parsing the header line or if there was no header line found in the file
	if err := vcfStreamer.ParseHeader("#CHROM"); err != nil {
		fmt.Printf("Encountered the following error while trying to parse the Header line of the vcf file being streamed in. Terminating program\n %s\n", err)
		os.Exit(1)
	} else if !vcfStreamer.Header_Found {
		fmt.Printf("Expected the input vcf file %s, to have a header line containing the string #CHROM. This line is essential to map the genotype calls to individuals. Please ensure that this line is in the file. Terminating program...\n", vcfStreamer.Filename)
		os.Exit(1)
	}

	// make a list of errors
	var err []error

	resultObj := Result{Errors: err, Samples: make(map[string]bool)}

	process_variant_stream(vcfStreamer, &resultObj)

	var error_encountered bool
	for _, msg := range resultObj.Errors {
		if msg != nil {
			fmt.Printf("Error: %s\n", msg)
			error_encountered = true
		}
	}
	if error_encountered {
		fmt.Println("Encountered the above errors while parsing through the vcf file stream. Terminating program...")
		os.Exit(1)
	}

	output_fh, open_err := os.Create(output_filepath)
	if open_err != nil {
		fmt.Printf("The following error was encountered while opening the file: %s", open_err)
	}

	buffered_writer := bufio.NewWriter(output_fh)

	writer(buffered_writer, resultObj)
}
