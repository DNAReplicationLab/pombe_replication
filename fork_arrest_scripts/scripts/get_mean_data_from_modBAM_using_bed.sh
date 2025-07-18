#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J getMeanDataModBamBed
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Given a mod bam file and a bed file, print an average brdu corresponding to each bed entry

# usage
#------
# sbatch get_mean_data_from_modBAM_using_bed.sh bed modbam
# bed and modbam refer to a bed file and a modified bam file respectively
# bed file should have at least 4 tab-separated columns with no column name: contig, start, end, read_id
# modbam file should be a modified bam file
# you can specify - as the bed file to read from standard input i.e. pipe bed file to standard input
# e.g. cat bed_file | bash get_mean_data_from_modBAM_using_bed.sh - modbam

# outputs
# -------
# sent to standard output

# stop execution if any command fails
set -e

# load packages
source load_package.sh -python

# check that the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
  echo >&2 "ERROR: incorrect number of arguments provided"
  echo >&2 "usage: bash get_mean_data_from_modBAM_using_bed.sh \$bed \$modbam"
  echo >&2 "bed and modbam refer to a bed file and a modified bam file respectively"
  echo >&2 "bed file should have at least 4 tab-separated columns with no column name: contig, start, end, read_id"
  echo >&2 "modbam file should be a modified bam file"
  echo >&2 "you can specify - as the bed file to read from standard input i.e. pipe bed file to standard input"
  echo >&2 "e.g. cat bed_file | bash get_mean_data_from_modBAM_using_bed.sh - modbam"
  exit 1
fi

# receive the input file
if [ "$1" == "-" ]; then

  # create a temporary directory
  temp_dir="${config[scratchDir]:-}"/tmp/"$(uuidgen)"
  mkdir -p "$temp_dir"

  # read from standard input and write to the temporary file
  input_file="$temp_dir"/input.bed
  cat > "$input_file"

else

  # check that the input file exist
  if [ ! -f "$1" ]; then
    echo >&2 "ERROR: $1 does not exist"
    exit 1
  fi

  input_file="$1"

fi

if [ ! -f "$2" ]; then
  echo >&2 "ERROR: $2 does not exist"
  exit 1
fi

if [ ! "$(< "$input_file" python validate_bed_format.py --allow-float-score --require-uuid)" == "valid" ]; then
  >&2 echo "Error: bed file format invalid!"
  exit 1;
fi

# NOTE: we make up a fork index just for the sake of the workflow below, it is not used anywhere.
# shellcheck disable=SC2016
< "$input_file" python print_valid_data_bed_lines.py |\
  awk 'BEGIN{OFS=" "}{print $1, $2, $3, $4, $4 "_" $1 "_" $2 "_" $3 "_unmapped_L_" $2 "_" $3}' |\
    sed '1icontig start end read_id alt_read_id' |\
      python get_raw_data_from_modBAM.py --piped-regions "$2" --alt-read-id-column |\
        sed '1idetectIndex\tposOnRef\tval' |\
          python get_mean_brdU_window.py --thres 0.5 |\
            awk -v IFS=" " -v OFS="\t" '{print $1, $2, $3, $4}' |\
              python split_fork_index_from_tsv.py --col 1 '$contig, $3, $4, $read_id, $2'

# delete the temporary directory if it was created
if [ -n "$temp_dir" ]; then
  rm -rf "$temp_dir"
fi