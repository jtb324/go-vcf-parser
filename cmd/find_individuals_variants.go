package cmd

import (
	"bufio"
	"fmt"
	internal "go-phers-parser/internal"
	"go-phers-parser/internal/files"
	"log/slog"
	"os"
	"slices"
	"strings"
	"time"
)

type SampleInfo struct {
	Score                 string
	PathogenicVariants    []string
	NonsynonymousVariants []string
	OtherVariants         []string
}

type SampleID struct {
	Index    int
	SampleID string
	Score    string
}

func read_samples_file(samples_filepath string, logger *slog.Logger) ([]string, []error) {
	logger.Info(fmt.Sprintf("Reading in all of the desired samples from the file %s\n", samples_filepath))
	var errors []error
	var samples []string

	samples_fh, samples_err := os.Open(samples_filepath)
	if samples_err != nil {
		errors = append(errors, fmt.Errorf("failed to open the file, %s. The following error was encountered, %s", samples_filepath, samples_err))
	} else {
		sample_scanner := bufio.NewScanner(samples_fh)
		for sample_scanner.Scan() {
			line := sample_scanner.Text()
			if strings.Contains(strings.ToLower(line), "grid") {
				// we can skip the header line if it exists
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
		logger.Info(fmt.Sprintf("Read %d samples in from the file, %s\n", len(samples), samples_filepath))
	}
	return samples, errors
}

func get_sample_col_indices(header_map map[string]int, samples []string, logger *slog.Logger) []SampleID {
	var sample_map []SampleID

	for sample_id, indx := range header_map {
		// Sometimes the id will have a score (either PheRS or case/control status) appended to the end. We can split the string to get this value
		split_id := strings.Split(sample_id, "_")

		// we can check if the sample id (split_id[0]) is in our samples array. If not then we skip that position
		if !slices.Contains(samples, split_id[0]) {
			continue
		}
		if len(split_id) == 2 {
			sample_map = append(sample_map, SampleID{Index: indx, SampleID: split_id[0], Score: split_id[1]})
		} else {
			sample_map = append(sample_map, SampleID{Index: indx, SampleID: split_id[0], Score: ""})
		}
	}
	logger.Info(fmt.Sprintf("Successfully mapped the indices for %d columns from the header", len(sample_map)))
	return sample_map
}

func check_for_alt_call(call string, reference_call_set map[string]bool) bool {
	_, ok := reference_call_set[call]
	// If the call is an alt call then the dictionary will return
	// false because it only contains keys for reference calls.
	// We can negate the false value so that it equals true
	// indicating that we did find an alternate allele
	return !ok
}

func find_col_indx(colname string, header_map map[string]int) (int, error) {
	col_indx, key_present := header_map[colname]

	if !key_present {
		return -1, fmt.Errorf("was not able to find the column, %s, in the header of the calls file. Please make sure that this column name was spelled exactly as it is found in the file", colname)
	}

	return col_indx, nil
}

func check_column_label(label string, values_of_interest []string) bool {
	var value_found bool
	for _, value := range values_of_interest {
		if strings.Contains(label, value) {
			value_found = true
			break
		}
	}
	return value_found
}

func initialize_sample_info(samples []SampleID) map[string]*SampleInfo {
	sampleInfo := make(map[string]*SampleInfo) // This will be our return value

	for _, obj := range samples {
		sampleInfo[obj.SampleID] = &SampleInfo{Score: obj.Score}
	}

	return sampleInfo
}

func parse_calls(calls_file string, samples []string, pathogenic_colname string, consequence_colname string, logger *slog.Logger) (map[string]*SampleInfo, []error) {
	var errors []error

	calls_fr := files.MakeFileReader(calls_file, 1024*1024)

	if calls_fr.Err != nil {
		fmt.Println(calls_fr.Err)
	}
	// lets defer the file closing
	// lets go ahead and parse through the calls_file to get the header
	err := calls_fr.ParseHeader("#CHROM")

	errors = append(errors, err)

	defer func() {
		for _, handle := range calls_fr.Handles {
			handle.Close()
		}
	}()

	// If we never found the header then we need to early exit. Other wise we will try to get an index that doesn't exist
	if !calls_fr.Header_Found {
		return nil, errors
	}
	// We need to find the columns for clinvar and the consequence columns

	clinVar_col_indx, clinvar_dict_err := find_col_indx(pathogenic_colname, calls_fr.Header_col_indx)

	consequence_col_indx, consequence_dict_err := find_col_indx(consequence_colname, calls_fr.Header_col_indx)

	if clinvar_dict_err != nil || consequence_dict_err != nil {
		errors = append(errors, clinvar_dict_err)
		errors = append(errors, consequence_dict_err)
		return nil, errors
	}
	// we also need to map the sample id columns
	sample_indices := get_sample_col_indices(calls_fr.Header_col_indx, samples, logger)

	sampleInfo := initialize_sample_info(sample_indices)

	// We also need to generate the set of reference calls so that we can compare our calls for that
	reference_call_strs := generate_reference_set()
	// This file has a header line so we first need to read in the indices for each column
	for calls_fr.FileScanner.Scan() {
		line := calls_fr.FileScanner.Text()
		// We assume the header line contains the phrase #CHROM because this is the output of the other program
		split_line := strings.Split(strings.TrimSpace(line), "\t")

		is_pathogenic := check_column_label(split_line[clinVar_col_indx], []string{"pathogenic", "likely_pathogenic"})
		is_nonsense_variant := check_column_label(split_line[consequence_col_indx], []string{"missense", "nonsynonymous"})

		for _, individual := range sample_indices {
			call := split_line[individual.Index]
			alternate_call := check_for_alt_call(call, reference_call_strs)
			// Now we can generate teh variant string that we are going to write to a file
			variantStr := fmt.Sprintf("%s:%s", split_line[2], call)
			individualInfo := sampleInfo[individual.SampleID]

			if is_pathogenic && alternate_call {
				individualInfo.PathogenicVariants = append(individualInfo.PathogenicVariants, variantStr)
			}

			if is_nonsense_variant && alternate_call {
				individualInfo.NonsynonymousVariants = append(individualInfo.NonsynonymousVariants, variantStr)
			}

			if !is_nonsense_variant && !is_pathogenic && alternate_call {
				individualInfo.OtherVariants = append(individualInfo.OtherVariants, variantStr)
			}

			// if check_for_alt_call(call, reference_call_strs) {
			// 	// We need to pull out the label for pathogenicity if that is present in the file
			// 	var pathogenic_label string
			// 	if pathogenic_label_present {
			// 		pathogenic_label = split_line[clinical_col_indx]
			// 	} else {
			// 		pathogenic_label = ""
			// 	}
			// 	// We need to keep track of the variants that an individual has. The variant is the first value in the split_line array.
			// 	if sampleStruct, ok := sampleInfo[individual.SampleID]; ok {
			// 		update_variant_status(sampleStruct, variantStr, pathogenic_label)
			// 	} else if check_pathogenic_label(pathogenic_label) {
			// 		variantList := []string{variantStr}
			// 		sampleInfo[individual.SampleID] = &SampleInfo{Score: individual.Score, PathogenicVariants: variantList, OtherVariants: []string{}}
			// 	} else {
			// 		nonPathogenicsList := []string{variantStr}
			// 		sampleInfo[individual.SampleID] = &SampleInfo{Score: individual.Score, PathogenicVariants: []string{}, OtherVariants: nonPathogenicsList}
			// 	}
			// }
		}
	}
	if calls_fr.FileScanner.Err() != nil {
		errors = append(errors, fmt.Errorf("encountered the following error while trying to scan through the calls file:  %s", calls_fr.FileScanner.Err()))
	}

	return sampleInfo, errors
}

func write_variants(writer *bufio.Writer, sample_variants map[string]*SampleInfo) {
	// lets build the header line

	header_str := "SAMPLE\tSCORE\tPATHOGENIC_VARIANTS\tNONSYNONYMOUS_VARIANTS\tOTHER_VARIANTS\n"

	writer.WriteString(header_str)

	sample_str := strings.Builder{}
	for sample_id, sampleInfoObj := range sample_variants {

		sample_str.WriteString(sample_id)

		pathogenicVarStr := strings.Join(sampleInfoObj.PathogenicVariants, ",")
		nonsynonymousVarStr := strings.Join(sampleInfoObj.NonsynonymousVariants, ",")
		otherVarStr := strings.Join(sampleInfoObj.OtherVariants, ",")

		// We can build the rest of the string appending the Score if there is one and the variants
		if sampleInfoObj.Score == "" {
			sample_str.WriteString(fmt.Sprintf("\t-\t%s\t%s\t%s", pathogenicVarStr, nonsynonymousVarStr, otherVarStr))
		} else {
			sample_str.WriteString(fmt.Sprintf("\t%s\t%s\t%s\t%s", sampleInfoObj.Score, pathogenicVarStr, nonsynonymousVarStr, otherVarStr))
		}
		sample_str.WriteString("\n")
	}

	writer.WriteString(sample_str.String())
	writer.Flush()
}

func FindSampleVariants(config internal.UserArgs, logger *slog.Logger) {
	start_time := time.Now()

	logger.Info(fmt.Sprintf("began the analysis at: %s\n", start_time.Format("2006-01-02@15:04:05")))
	// read in the appropriate CLI flags

	var samples []string
	var sample_file_err []error
	if config.PhenoFilePath == "" {
		logger.Error("No file contained the list of cases was provided. Please make sure you provide a file where the first column list all of the cases in the network to pul variants for")
	} else if config.PhenoFilePath != "" {
		// process the samples file
		samples, sample_file_err = read_samples_file(config.PhenoFilePath, logger)
		if sample_file_err != nil {
			logger.Error(fmt.Sprintf("Encountered the following errors while trying to read in samples from the file %s\n", config.PhenoFilePath))
			for msg_indx, msg := range sample_file_err {
				logger.Error(fmt.Sprintf("Error Msg %d:\n %s", msg_indx, msg))
			}
			os.Exit(1)
		}
	}
	// now we can parse through the output file for variants of interest

	// Create the scanner to read the calls file with a custom buffer

	sample_variants, errs := parse_calls(config.CallsFile, samples, config.ClinvarColumnName, config.ConsequenceCol, logger)

	var parsing_err_encountered bool
	for _, err_msg := range errs {
		if err_msg != nil {
			logger.Error(fmt.Sprintf("Error Msg:\n%s\n", err_msg))
			parsing_err_encountered = true
		}
	}
	if parsing_err_encountered {
		logger.Info("Terminating program because of the above errors...")
		os.Exit(1)
	}

	logger.Info(fmt.Sprintf("Identified variants for %d samples", len(sample_variants)))

	output_fh, output_err := os.Create(config.OutputFilepath)

	if output_err != nil {
		logger.Error(fmt.Sprintf("Encountered the following error while trying to open the output file, %s.\n %s", config.OutputFilepath, output_err))
		os.Exit(1)
	}

	defer output_fh.Close()

	writer := bufio.NewWriter(output_fh)
	logger.Info(fmt.Sprintf("Writing output to the file: %s", config.OutputFilepath))
	write_variants(writer, sample_variants)

	end_time := time.Now()

	logger.Info(fmt.Sprintf("finished analysis at: %s", end_time.Format("2006-01-02@15:04:05")))

	duration := end_time.Sub(start_time)

	logger.Info(fmt.Sprintf("total analysis time: %s", duration.String()))
}
