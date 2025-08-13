package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var RootCmd = &cobra.Command{
	Use:   "go-variant-parser-{platform}",
	Short: "CLI tool to read through genetic sequencing data and add information about who are the carriers and what are their phers values",
	// Run: func(cmd *cobra.Command, args []string) {
	// },
	// PersistentPreRun: func(cmd *cobra.Command, args []string) {
	// 	f, err := os.Create("test_cpu_profiling.prof")
	// 	if err != nil {
	// 		fmt.Fprintf(os.Stderr, "Error creating CPU profile: %v\n", err)
	// 		os.Exit(1)
	// 	}
	// 	if err := pprof.StartCPUProfile(f); err != nil {
	// 		fmt.Fprintf(os.Stderr, "Error starting CPU profile: %v\n", err)
	// 		os.Exit(1)
	// 	}
	//
	// 	// If the command is interrupted before the end (ctrl-c), flush the
	// 	// profiling files
	// 	c := make(chan os.Signal, 1)
	// 	signal.Notify(c, os.Interrupt)
	// 	go func() {
	// 		<-c
	// 		f.Close()
	// 		os.Exit(0)
	// 	}()
	// },
	// PersistentPostRun: func(cmd *cobra.Command, args []string) {
	// 	pprof.StopCPUProfile()
	// },
}

func init() {
	RootCmd.CompletionOptions.DisableDefaultCmd = true
}

func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Fprintf(os.Stderr, "Encountered the following error while instantating the application '%s'\n", err)
		os.Exit(1)
	}
}
