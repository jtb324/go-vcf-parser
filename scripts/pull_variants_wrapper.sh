#!/bin/bash
# We are going to make some of the font bold so it stands out. We can use tput so that this works on more systems than just the stand escape sequences
bold=$(tput bold)
normal=$(tput sgr0)
underline=$(tput smul)
no_underline=$(tput rmul)

if [ "$#" -ne 7 ]; then
    echo "${bold}ERROR:${normal} No arguments provided to the run script"
    echo "Expected command to be 'bash pull_polg_variants.sh {chrX:start-end} {samples_filepath} {sequencing_file} {annotation_filepath} {phers_samples} {annotation_cols_to_keep} {output_path}"
else
    REGION=$1
    SAMPLES_FILE=$2
    SEQUENCING_FILE=$3
    ANNO_FILE=$4
    PHERS_FILE=$5
    ANNO_COLS=$6 # This argument should be a comma separated list. Ex: col1,col2,col3. And the columns need to be spelled the exact same way as they are in the annotation file
    OUTPUT_PATH=$7
    # Lets print out information to the log file
    printf '=%.0s' {1..40}
    echo -e "\n"
    echo "The following arguments were provided to the script:"
    echo "${bold}REGION${normal}=${REGION}"
    echo "${bold}SAMPLES_FILE${normal}=${SAMPLES_FILE}"
    echo "${bold}SEQUENCING_FILE${normal}=${SEQUENCING_FILE}"
    echo "${bold}ANNO_FILE${normal}=${ANNO_FILE}"
    echo "${bold}PHERS SCORES${normal}=${PHERS_FILE}"
    echo "${bold}ANNOTATION COLS TO KEEP${normal}=${ANNO_COLS}" 
    echo "${bold}OUTPUT PATH${normal}=${OUTPUT_PATH}"
    echo -e "\n"
    printf '=%.0s' {1..40}
    echo -e "\n"
    
    echo "${underline}command being run:${no_underline}"
    echo "bcftools view -r ${REGION} -S ${SAMPLES_FILE} -Ov ${SEQUENCING_FILE} | go-variant-parser-linux-amd64 pull-variants -a ${ANNO_FILE} -c ${ANNO_COLS} -o ${OUTPUT_PATH} --maf-threshold 0.1 -s ${PHERS_FILE} -r ${REGION}"
    echo -e "\n"

    # Now we can actually run the command we want to
    bcftools view -r "$REGION" -S "$SAMPLES_FILE" -Ov "$SEQUENCING_FILE" | go-variant-parser-linux-amd64 pull-variants -a "$ANNO_FILE" -c "$ANNO_COLS" -o "$OUTPUT_PATH" --maf-threshold 0.1 -s "$PHERS_FILE" -r "$REGION"
fi

