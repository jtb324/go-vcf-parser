package cmd

import (
	"bufio"
	"fmt"
	"go-phers-parser/internal/files"
	"maps"
	"os"
	"slices"
	"strings"

	"github.com/spf13/cobra"
)

var variant_calls = &cobra.Command{
	Use:   "find_variant_calls",
	Short: "find the individuals with variant calls for a site of interest. Expects vcf input to be streamed in from bcftools",
	Run:   find_variant_calls,
}

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
	if count, value_exist := genotype_counts[call]; value_exist {
		count++
	} else {
		genotype_counts[call] = 1
	}
}

func process_line(line string, header_map map[int]string, resultsObj Result) {
	variantCallsObj := VariantCalls{VariantCarriers: make(map[string]string), GenotypeCounts: make(map[string]int)}

	split_line := strings.Split(strings.TrimSpace(line), "\t")

	// We can add the variant string here
	variantCallsObj.VariantInfo = split_line[0:3]

	// We will need to generate the reference calls for comparison
	ref_call_set := generate_reference_set()
	// We can iterate over each call
	for indx, calls := range split_line[9:] {
		if check_alt_call(calls, ref_call_set) {
			// We can add the id and the call to the carriers map
			id := header_map[indx]
			variantCallsObj.VariantCarriers[id] = calls
			// Then we can also save the carrier ids we found. We will use
			// this list to create the header for the output file later
			resultsObj.Samples[id] = true // This is how you use a set in Go. Its the same as a map
		}
		update_genotype_count(calls, variantCallsObj.GenotypeCounts)
	}
	fmt.Printf("Identified %d individuals who were either heterozygous or homozygous alt for the variant %s", len(variantCallsObj.VariantCarriers), variantCallsObj.VariantInfo[2])

	resultsObj.Variants = append(resultsObj.Variants, variantCallsObj)
}

func writer(writer *bufio.Writer, results Result) {
	// get a list of all the samples we need to put in the header
	sample_list := results.generate_sample_list()
	// Create the header string
	header_str := strings.Builder{}
	header_str.WriteString("CHROM\tPOS\tID\t")
	header_str.WriteString(fmt.Sprintf("%s\t", strings.Join(sample_list, "\t")))

	writer.WriteString(header_str.String())
}

// This function is used to find all the individuals with variant calls for a site of interest.
// It expects to have input streamed in from bcftools
func find_variant_calls(cmd *cobra.Command, args []string) {
	output_filepath, _ := cmd.Flags().GetString("output")
	// we need to create the reader
	vcfStreamer := files.NewVcfStreamer(5012 * 5012)

	if err := vcfStreamer.Initialize(); err != nil {
		fmt.Printf("Encountered the following error while reading through the streamed vcf file header")
	}

	// make a list of errors
	var err []error

	resultObj := Result{Errors: err, Samples: make(map[string]bool)}

	for vcfStreamer.Next_line != "" {
		process_line(vcfStreamer.Next_line, vcfStreamer.Header_col_mappings, resultObj)
		vcfStreamer.ReadNextLine()
	}

	if vcfStreamer.CheckErrs() != nil {
		fmt.Println("There was an error")
	}

	output_fh, open_err := os.Create(output_filepath)
	if open_err != nil {
		fmt.Printf("The following error was encountered while opening the file: %s", open_err)
	}

	buffered_writer := bufio.NewWriter(output_fh)

	writer(buffered_writer, resultObj)
}

func init() {
	RootCmd.AddCommand(variant_calls)
	variant_calls.Flags().StringP("output", "o", "test_output.txt", "Filepath to write the output file to.")
	variant_calls.Flags().Int64P("buffersize", "b", 5012*5012, "buffersize to use while reading through the streamed input data. Default: 5012*5012 bytes")
}
