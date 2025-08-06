BINARY_NAME=go-variant-parser
BUILD_DIR=./build

.PHONY: help
help:
		@echo "Usage:"
		@sed -n 's/^##//p' ${MAKEFILE_LIST} | column -t -s ':' | sed -e 's/^/ /'

## build: Build the application into a folder called "build". Binaries will be made for linux-amd64, darwin-arm64, and darwin-amd64
.PHONY: build
build:
		@mkdir -p ${BUILD_DIR}
		@echo "creating binaries for linux-amd64, darwin-amd64, and darwin-arm64"
		GOARCH=amd64 GOOS=darwin go build -o ${BUILD_DIR}/${BINARY_NAME}-darwin-amd64 .
		GOARCH=amd64 GOOS=linux go build -o ${BUILD_DIR}/${BINARY_NAME}-linux-amd64 .
		GOARCH=arm64 GOOS=darwin go build -o ${BUILD_DIR}/${BINARY_NAME}-darwin-arm64 .

.PHONY: confirm
confirm: 
		@echo -n 'Please confirm that you wish to remove the build directory, ${BUILD_DIR}. [y/N] ' && read ans && [ $${ans:-N} = y ] # We use the ans value but if the ans is empty then we will use N as a default

## clean: clean up the build folder
.PHONY: clean
clean: confirm
		@echo "removing the directory ${BUILD_DIR}"
		@rm -rf "${BUILD_DIR}"
