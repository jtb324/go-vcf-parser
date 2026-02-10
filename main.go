package main

import (
	"context"
	"log"
	"os"

	cmd_commands "go-phers-parser/cmd"

	"github.com/urfave/cli/v3"
)

func main() {
	cmd := &cli.Command{
		Name:  "go-vcf-parser",
		Usage: "A small go utility to parse vcf files",
		Commands: []*cli.Command{
			{
				Name:  "pull-variants",
				Usage: "pull variants for the specified region",
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:    "anno-file",
						Aliases: []string{"a"},
						Usage:   "Filepath to an annotation file (currently on supports VEP so that there is a canocial colum that we can use to avoid duplicates and only look at the cannocial transcript).",
					},
					&cli.StringFlag{
						Name:    "keep-cols",
						Aliases: []string{"c"},
						Usage:   "Columns in the annotation file to keep while it is being read in.",
					},
					&cli.StringFlag{
						Name:    "samples",
						Aliases: []string{"s"},
						Usage:   "File where the first column contains the samples that we wish to keep in the analysis",
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
					&cli.StringFlag{
						Name:  "log-filepath",
						Value: "test.log",
						Usage: "Filepath to write the log file to.",
					},
					&cli.FloatFlag{
						Name:  "maf-threshold",
						Value: 0.1,
						Usage: "Minor allele frequency cap to filter output so that only variants below this threshold are returned",
					},
					&cli.IntFlag{
						Name:    "buffersize",
						Aliases: []string{"b"},
						Value:   5012 * 5012,
						Usage:   "buffersize to use while reading through the streamed input data. Default: 5012&**2 bytes",
					},
				},
				Action: func(ctx context.Context, cmd *cli.Command) error {
					pull_vars_args := cmd_commands.PullVariantsArgs{
						AnnoFile:        cmd.String("anno-file"),
						ColsToKeep:      cmd.String("keep-cols"),
						SamplesFilepath: cmd.String("samples"),
						OutputFile:      cmd.String("output"),
						LogFilePath:     cmd.String("log-filepath"),
						MafCap:          cmd.Float("maf_threshold"),
						Region:          cmd.String("region"),
						Buffersize:      cmd.Int("buffersize"),
					}

					cmd_commands.PullVariants(pull_vars_args)

					return nil
				},
			},
			{
				Name:  "find-all-carriers",
				Usage: "find the individuals with variant calls for a site of interest. Expects vcf input to be streamed in from bcftools",
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:    "output",
						Aliases: []string{"o"},
						Usage:   "Filepath ot write the output file to.",
					},
					&cli.IntFlag{
						Name:    "buffersize",
						Aliases: []string{"b"},
						Value:   5012 * 5012,
						Usage:   "buffersize to use while reading through the streamed input data. Default: 5012&**2 bytes",
					},
					&cli.StringFlag{
						Name:  "sample-exclusion-string",
						Usage: "List of comma-separated substrings that may indicate if a sample should be excluded from the analysis. This situation can arise if the reference panel controls were kept in the vcf or if invalid samples are present. This code can filter out those individuals by seeing if the substring is present in the ID. This list should not have spaces between the strings",
					},
				},
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
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:    "calls-file",
						Aliases: []string{"c"},
						Usage:   "output tab separated file from the pull-variants command. Each row should be a variant and each column willl be an individual with tthe phers score (if the score was provided in that command)",
					},
					&cli.StringFlag{
						Name:    "samples-list",
						Aliases: []string{"s"},
						Usage:   "list of sample ids to find all the variants for. This list shoudl be comma separated with spaces no spaces between the ids",
					},
					&cli.StringFlag{
						Name:    "samples-file",
						Aliases: []string{"S"},
						Usage:   "filepath to a tab separated text file that has the ssampels we wish to keep. The first column is expected to be a list of grids and the file amy or may not have a header",
					},
					&cli.StringFlag{
						Name:    "output",
						Aliases: []string{"o"},
						Value:   "test_output.txt",
						Usage:   "Filepath to write the output file to.",
					},
					&cli.StringFlag{
						Name:  "clinvar-col",
						Usage: "column label of the clinical annotations column. These annotations can come fro VEP or manual annotations.",
					},
					&cli.StringFlag{
						Name:  "consequence-col",
						Usage: "column label of the consequences columns. This column shoudl contain values like 'intron_variant' or 'missense_variant', etc...",
					},
					&cli.StringFlag{
						Name:  "log-filepath",
						Value: "test.log",
						Usage: "Filepath to write the log messages to",
					},
				},
				Action: func(ctx context.Context, cmd *cli.Command) error {
					userArgs := cmd_commands.RuntimeConfig{
						CallsFile:         cmd.String("calls-file"),
						SamplesList:       cmd.String("samples-list"),
						SamplesFilepath:   cmd.String("samples-file"),
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
		},
	}
	if err := cmd.Run(context.Background(), os.Args); err != nil {
		log.Fatal(err)
	}
}
