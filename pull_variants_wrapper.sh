#!/bin/bash
# We are going to make some of the font bold so it stands out. We can use tput so that this works on more systems than just the stand escape sequences


usage() {
  echo "Usage: {scriptname} [-r] [-s] [-i] [-a] [-p] [-k] [-m] [-o] [-c] [-C]"
  printf "=%.0b" {1..70}
  printf "\n"
  echo "Options:"
  echo "  -h    show the help message"
  echo "Required Flags:"
  echo "  -r    region of the genome to search for"
  echo "  -s    tab separated text file containing all the samples that we want to pull the sequencing calls for. This file should be formatted for bcftools -S flag"
  echo "  -i    sequencing vcf file to use as input. This will be passed to bcftools"
  echo "  -a    annotation file from VEP that will be search through"
  echo "  -p    tab separated text file where the first column is the individuals IDs and the second column is a score or case/control status. Expected to not have a header"
  echo "  -k    Columns that the program will keep when it parses through the annotation file"
  echo "  -m    minor allele frequency cap that defines the maximum minor allele frequency that we will use to parse the sequencing file."
  echo "  -c    column index (0-based) of the clinical annotations within the annotation file"
  echo "  -n    a file of cases in the network of interest to pull the sequencing for. This file should have 1 column and may or may not have a header. This file should be a subset of the '-s' samples file"
  echo "  -o    Output filepath with a prefix to write to. This argument should not end with any file extension and it should end in a directory path. The program will append the suffixes  *_all_network_id_variants.txt and *_cases_in_network_variants.txt to create 2 files"
  echo "  -v    name of the column that has variant conseqeunces such as 'missense' or 'intro_variant'"
}
optstring="hr:s:i:a:p:k:o:m:c:n:v:"

# while getopts ":h:r:s:i:a:p:k:o:m:c:n:v" opt; do
while getopts ${optstring} opt; do
  case ${opt} in
  h) 
    usage
    exit 0
    ;;
  r) REGION=${OPTARG} ;;
  s) SAMPLES_FILE=${OPTARG} ;;
  i) SEQUENCING_FILE=${OPTARG} ;;
  a) ANNO_FILE=${OPTARG} ;;
  p) PHERS_FILE=${OPTARG} ;;
  k) ANNO_COLS=${OPTARG} ;;
  o) OUTPUT_PATH=${OPTARG} ;;
  m) MAF_THRESHOLD=${OPTARG} ;;
  c) CLINVAR_COL_INDX=${OPTARG} ;;
  n) NETWORK_CASE_FILE=${OPTARG} ;;
  v) CONSEQUENCE_COL_STRING=${OPTARG} ;;
  \?)
    echo "Invalid option: -$OPTARG" 1>&2
    echo "Usage Menu shown below:"
    printf "\n"
    usage
    exit 1
    ;;
  : ) 
    echo "Invalid option: $OPTARG requires an argument" 1>&2
    usage
    exit 1
    ;;
  esac
done

# We can first check if the bcftools command is avaliable
if ! command -v bcftools &>/dev/null || ! command -v go-variant-parser-linux-amd64 &>/dev/null; then
  echo "Unable to find the program bcftools or the executable go-variant-parser-linux-amd64. Please make sure that these programs are installed and placed in your PATH so that they can be detected"
  exit 1
fi
# If a MAF_THRESHOLD wasn't passed then we are going to default to 0.1
if [[ $MAF_THRESHOLD ]]; then MAF=$MAF_THRESHOLD; else MAF=0.1; fi
# If it is present then we can continue on with the script
# Lets print out information to the log file
printf '=%.0s' {1..60}
printf "\n"
echo "The following arguments were provided to the script:"
echo "REGION=${REGION}"
echo "SAMPLES_FILE=${SAMPLES_FILE}"
echo "SEQUENCING_FILE=${SEQUENCING_FILE}"
echo "ANNO_FILE=${ANNO_FILE}"
echo "PHERS SCORES=${PHERS_FILE}"
echo "ANNOTATION COLS TO KEEP=${ANNO_COLS}"
echo "OUTPUT PREFIX PATH=${OUTPUT_PATH}"
echo "MAXMIUM ALLELE FRQEUENCY THRESHOLD=${MAF}"
echo "CLINVAR LABEL INDEX=${CLINVAR_COL_INDX}"
echo "NETWORK CASES FILE=${NETWORK_CASE_FILE}"
printf '=%.0s' {1..60}
printf "\n%.0s" {1..5}

## We need to make the output path for both the first go command and the second
OUTPUT_GO_CMD1=${OUTPUT_PATH}_all_network_id_variants.txt
OUTPUT_GO_CMD2=${OUTPUT_PATH}_cases_in_network_variants.txt

echo "Reading in annotations for the region ${REGION} and pulling variants for the samples in the samples file, ${SAMPLES_FILE}"
printf "\n"
printf '=%.0s' {1..60}
printf "\n"
#
echo "${underline}command being run:${no_underline}"
echo "bcftools view -r ${REGION} -S ${SAMPLES_FILE} -Ov ${SEQUENCING_FILE} | go-variant-parser-linux-amd64 pull-variants -a ${ANNO_FILE} -c ${ANNO_COLS} -o ${OUTPUT_GO_CMD1} --maf-threshold 0.1 -s ${PHERS_FILE} -r ${REGION}"
printf "\n%.0s" {1..3}
# #
# # # Now we can actually run the command we want to
bcftools view -r "$REGION" -S "$SAMPLES_FILE" -Ov "$SEQUENCING_FILE" | go-variant-parser-linux-amd64 pull-variants -a "$ANNO_FILE" -c "$ANNO_COLS" -o "$OUTPUT_GO_CMD1" --maf-threshold "$MAF" -s "$PHERS_FILE" -r "$REGION"
# #
printf "=%.0s" {1..60}
printf "\n"
echo "Finding all the variants that the samples in the file, ${NETWORK_CASE_FILE}, have using the output file ${OUTPUT_GO_CMD1}"
printf "\n"
printf '=%.0s' {1..60}
printf "\n"
# Now we can find the other variants that the cases of the network have. We
# have to use the output of the first go command as input
go-variant-parser-linux-amd64 view-sample-variants -c "$OUTPUT_GO_CMD1" --clinvar-label-present --clinvar-col "$CLINVAR_COL_INDX" -S "$NETWORK_CASE_FILE" -o "$OUTPUT_GO_CMD2"
