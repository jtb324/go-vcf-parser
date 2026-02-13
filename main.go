package main

import (
	"context"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"

	cmd_commands "go-phers-parser/cmd"
	"go-phers-parser/internal"

	"github.com/urfave/cli/v3"
)

func main() {
	// we are going to define our flag arrays here
	pull_var_flags := []cli.Flag{
		&cli.StringFlag{
			Name:    "anno-file",
			Aliases: []string{"a"},
			Usage:   "Filepath to an annotation file (currently on supports VEP so that there is a canocial colum that we can use to avoid duplicates and only look at the cannocial transcript).",
		},
		&cli.StringFlag{
			Name:    "pheno-file",
			Aliases: []string{"p"},
			Usage:   "Filepath to a tab separated file where the first column are ids and the second column is the case/control status. This file can have a header with the columns 'GRID' and 'Status' or it can have no header",
		},
		&cli.StringFlag{
			Name:    "keep-cols",
			Aliases: []string{"c"},
			Usage:   "Columns in the annotation file to keep while it is being read in.",
		},
		&cli.StringFlag{
			Name:    "output",
			Aliases: []string{"o"},
			Value:   "test_output.txt",
			Usage:   "Filepath to write the output file to.",
		},
		&cli.StringFlag{
			Name:    "region",
			Aliases: []string{"r"},
			Usage:   "region of the chromosome that we are intereseted in search for. This regions should have the form chrX:start-end. We will use this region to filter which annotations we wish to save in memory",
		},

		&cli.FloatFlag{
			Name:  "maf-threshold",
			Value: 0.1,
			Usage: "Minor allele frequency cap to filter output so that only variants below this threshold are returned",
		},
	}

	find_all_carriers_flags := []cli.Flag{
		&cli.StringFlag{
			Name:  "sample-exclusion-string",
			Usage: "List of comma-separated substrings that may indicate if a sample should be excluded from the analysis. This situation can arise if the reference panel controls were kept in the vcf or if invalid samples are present. This code can filter out those individuals by seeing if the substring is present in the ID. This list should not have spaces between the strings",
		},
	}

	pull_sample_variants := []cli.Flag{
		&cli.StringFlag{
			Name:  "clinvar-col",
			Usage: "column label of the clinical annotations column. These annotations can come fro VEP or manual annotations.",
		},
		&cli.StringFlag{
			Name:  "consequence-col",
			Usage: "column label of the consequences columns. This column shoudl contain values like 'intron_variant' or 'missense_variant', etc...",
		},
	}

	cmd := &cli.Command{
		Name:  "go-vcf-parser",
		Usage: "A small go utility to parse vcf files",
		// define global flags for all commands
		Flags: []cli.Flag{
			&cli.IntFlag{
				Name:    "buffersize",
				Aliases: []string{"b"},
				Value:   5012 * 5012,
				Usage:   "buffersize to use while reading through the streamed input data. Default: 5012**2 bytes",
			},
			&cli.StringFlag{
				Name:  "log-filepath",
				Value: "test.log",
				Usage: "Filepath to write the log file to.",
			},
			&cli.StringFlag{
				Name:    "output",
				Aliases: []string{"o"},
				Value:   "test_output.txt",
				Usage:   "Filepath to write the output file to. If running subcommands individually then this should be a full file path with a suffix. If you are running the pipeline command then this value should only be the output prefix.",
			},
		},
		Commands: []*cli.Command{
			{
				Name:  "pull-variants",
				Usage: "pull variants for the specified region",
				Flags: pull_var_flags,
				Action: func(ctx context.Context, cmd *cli.Command) error {
					pull_vars_args := internal.UserArgs{
						AnnoFile:      cmd.String("anno-file"),
						ColsToKeep:    cmd.String("keep-cols"),
						PhenoFilePath: cmd.String("pheno-file"),
						OutputFile:    cmd.String("output"),
						LogFilePath:   cmd.String("log-filepath"),
						MafCap:        cmd.Float("maf_threshold"),
						Region:        cmd.String("region"),
						Buffersize:    cmd.Int("buffersize"),
					}

					cmd_commands.PullVariants(pull_vars_args)

					return nil
				},
			},
			{
				Name:  "find-all-carriers",
				Usage: "find the individuals with variant calls for a site of interest. Expects vcf input to be streamed in from bcftools",
				Flags: find_all_carriers_flags,
				Action: func(ctx context.Context, cmd *cli.Command) error {
					output_path := cmd.String("output")
					buffersize := cmd.Int("buffersize")
					sample_exclusion := cmd.String("sample-exclusion-string")

					cmd_commands.FindAllCarrierCalls(output_path, buffersize, sample_exclusion)

					//TODO: Need to update the FindAllCarrierCalls to return an error
					return nil
				},
			},
			{
				Name:  "view-sample-variants",
				Usage: "grab the variants that samples of interest have. This command uses the output from the pull-variants command",
				Flags: pull_sample_variants,
				Action: func(ctx context.Context, cmd *cli.Command) error {
					userArgs := internal.UserArgs{
						CallsFile:         cmd.String("calls-file"),
						PhenoFilePath:     cmd.String("pheno-file"),
						OutputFilepath:    cmd.String("output"),
						ClinvarColumnName: cmd.String("clinvar-col"),
						ConsequenceCol:    cmd.String("consequence-col"),
						LogfilePath:       cmd.String("log-filepath"),
					}

					cmd_commands.FindSampleVariants(userArgs)

					//TODO: Need to update the FindSampleVariants to return an error
					return nil
				},
			},
			{
				Name:  "run-pipeline",
				Usage: "This subcommand serves as a pipeline that connects the pull-variants subcommand with the view-sample-variants subcommand. So that users can run both together if they wish to. To run this we are assuming that the input sequencing file is being piped through bcftools",
				// Now we can appened the subcommand flags to this pipeline
				Flags: append(append([]cli.Flag{}, pull_var_flags...), pull_sample_variants...),
				Action: func(ctx context.Context, cmd *cli.Command) error {

					start_time := time.Now()

					fmt.Printf("began the analysis at: %s\n", start_time.Format("2006-01-02@15:04:05"))

					// We need to first make sure that the output file has no suffix (meaning it is only a prefix)
					userProvidedOutput := cmd.String("output")
					final_output_prefix := strings.TrimSuffix(userProvidedOutput, filepath.Ext(userProvidedOutput))

					output_file1 := fmt.Sprintf("%s_all_network_id_variants.txt", final_output_prefix)

					output_file2 := fmt.Sprintf("%s_cases_in_network_variants.txt", final_output_prefix)

					userArgs := internal.UserArgs{
						AnnoFile:          cmd.String("anno-file"),
						ColsToKeep:        cmd.String("keep-cols"),
						OutputFile:        output_file1,
						LogFilePath:       cmd.String("log-filepath"),
						MafCap:            cmd.Float("maf_threshold"),
						Region:            cmd.String("region"),
						Buffersize:        cmd.Int("buffersize"),
						CallsFile:         output_file1,
						PhenoFilePath:     cmd.String("pheno-file"),
						OutputFilepath:    output_file1,
						ClinvarColumnName: cmd.String("clinvar-col"),
						ConsequenceCol:    cmd.String("consequence-col"),
						LogfilePath:       cmd.String("log-filepath"),
					}

					fmt.Printf("Reading in annotations for the region %s and pulling variants for the samples in the samples file, %s\n", userArgs.Region, userArgs.PhenoFilePath)

					cmd_commands.PullVariants(userArgs)

					//lest make sure that the output file is right now
					userArgs.OutputFile = output_file2

					cmd_commands.FindSampleVariants(userArgs)

					end_time := time.Now()

					fmt.Printf("finished analysis at: %s\n", end_time.Format("2006-01-02@15:04:05"))

					duration := end_time.Sub(start_time)

					fmt.Printf("total analysis time: %s\n", duration.String())
					//TODO: Need to update the FindSampleVariants to return an error
					return nil
				},
			},
		},
	}
	if err := cmd.Run(context.Background(), os.Args); err != nil {
		log.Fatal(err)
	}
}
