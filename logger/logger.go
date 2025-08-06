package logger

import (
	"fmt"
	"log/slog"
	"os"
)

// In the future we are going to try to implemement a more verbose logging level not just INFO or DEBUG
const (
	LevelVerbose = slog.Level(-2)
)

func CreateLogger(loglevel string, logFilePath string) (*slog.Logger, *slog.Logger) {
	// we can set the log level based on user input
	curr_log_level := &slog.LevelVar{}

	switch loglevel {
	case "debug":
		curr_log_level.Set(slog.LevelDebug)
	case "info":
		curr_log_level.Set(slog.LevelInfo)
	default:
		fmt.Fprintf(os.Stderr, "Did not recognize the logging level of %s", loglevel)

	}

	opts := &slog.HandlerOptions{
		AddSource: true,
		Level:     curr_log_level,
	}

	consoleLogger := slog.New(slog.NewTextHandler(os.Stdout, opts))

	var fileLogger *slog.Logger

	if logFilePath != "" {
		log_fh, file_err := os.Open(logFilePath)

		if file_err != nil {
			consoleLogger.Error("Unable to open the logging file. Logging will not be written to a file only STDOUT")
			fileLogger = nil
		} else {
			fileLogger = slog.New(slog.NewTextHandler(log_fh, opts))
		}
	} else {
		consoleLogger.Info("No filepath provided for the logging file. Only writting loig messages to STDOUT")
	}

	return consoleLogger, fileLogger
}
