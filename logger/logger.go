package logger

import (
	"log/slog"
	"os"
)

// In the future we are going to try to implemement a more verbose logging level not just INFO or DEBUG
const (
	LevelVerbose = slog.Level(-2)
)

func CreateLogger(loglevel int, logFilePath string) *slog.Logger {
	// we can set the log level based on user input
	curr_log_level := &slog.LevelVar{}

	switch loglevel {
	case 0:
		curr_log_level.Set(slog.LevelDebug)
	case 1:
		curr_log_level.Set(slog.LevelInfo)
	default:
		curr_log_level.Set(LevelVerbose)
	}

	opts := &slog.HandlerOptions{
		AddSource: true,
		Level:     curr_log_level,
	}

	consoleLogger := slog.New(slog.NewTextHandler(os.Stdout, opts))

	return consoleLogger
}
