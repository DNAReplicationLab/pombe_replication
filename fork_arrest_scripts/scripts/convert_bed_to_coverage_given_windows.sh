#!/bin/bash

# goal
# ====
# Given a bed file, calculate the coverage of the bed file at each given window.
# We require 100% coverage of the window by a bed entry or 100% coverage of the bed entry by the window to count it.

# usage
# =====
# bash convert_bed_to_coverage_given_windows.sh [-i] window_file bed_file fasta_index_file > output_file.bedgraph
# Do not use slurm here as the script is very fast.
#     window_file: a bed file of windows, must have only 3 columns.
#     bed_file: the bed file which is intersected with the windows to calculate the coverage.
#     fasta_index_file: the fasta index file for the genome.
#    -i: (invert) instead of requiring 100% coverage of the window by a bed entry, require 100% coverage of the bed entry by the window.
#    -i is optional.

# logic
# =====
# Given a window i and a bed entry j,
# if 100% coverage of i by j or 100% coverage of j by i if -i is used, increment the coverage of i by 1.

# output
# ======
# The output is to stdout, is space-delimited with no headers and four columns:
# contig, start, end, coverage
# each contig, start, end entry is from the window_file.
# Output can be redirected to a file.
# No header is printed.

# fail if any command fails
set -e

# parse options
invert="-f"
while getopts "i" opt; do
  case $opt in
    i)
      invert="-F"
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# shift the options off the arguments
shift $((OPTIND-1))

# check that the correct number of arguments were provided
if [ $# -lt 3 ]
then
    >&2 echo "Error: not enough arguments were provided"
    >&2 echo "Goal: Given a bed file, calculate the coverage of the bed file at each given window."
    >&2 echo "Usage: bash convert_bed_to_coverage.sh [-i] window_file bed_file fasta_index_file > output_file.bedgraph"
    >&2 echo "    window_file: a bed file of windows with only 3 columns."
    >&2 echo "    bed_file: the bed file which is intersected with the windows to calculate the coverage."
    >&2 echo "    fasta_index_file: the fasta index file for the genome (can be generated with samtools faidx). "
    >&2 echo "   -i: (optional, invert) instead of requiring 100% coverage of the window by a bed entry, require 100% coverage of the bed entry by the window."
    exit 1;
fi

# get the arguments
window_file=$1
bed_file=$2
fasta_index_file=$3

# check that the input files exist
if [ ! -f "$bed_file" ] || [ ! -f "$window_file" ] || [ ! -f "$fasta_index_file" ]
then
    >&2 echo "Error: the input files do not exist"
    exit 1;
fi

# load configuration variables, git labels, and bedtools
source load_package.sh -bedtools -python
source load_git_repo_labels.sh
source config.sh

# both bed files must have only three columns
if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: bed file is invalid."
  exit 1;
fi
if [ ! "$(< "$window_file" python validate_bed_format.py --max-three-columns)" == "valid" ]; then
  >&2 echo "Error: window file must have exactly three valid columns."
  exit 1;
fi

# both bed files must have contigs only in the fasta index file
if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fasta_index_file")" == "valid" ]; then
  >&2 echo "Error: bed file must have contigs only in the fasta index file."
  exit 1;
fi

if [ ! "$(< "$window_file" python validate_bed_against_fai.py "$fasta_index_file")" == "valid" ]; then
  >&2 echo "Error: bed file must have contigs only in the fasta index file."
  exit 1;
fi

# sort the bed file by contig and position, and retain only the first three columns
temp_file=$(mktemp)
bedtools sort -i "$bed_file" -faidx "$fasta_index_file" | awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1, $2, $3}}'  > "$temp_file"

# print the script information
echo "# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>";
echo "# script: $0";
echo "# arguments: $*";

# calculate the coverage of the windows by the forks and print four columns: contig, start, end, coverage
# if the invert flag is set, calculate the coverage of the bed file by the windows
bedtools coverage -a "$window_file" -b "$temp_file" -sorted -g "$fasta_index_file" $invert 1.0 |\
    awk 'BEGIN{OFS=" "}{print $1, $2, $3, $4}'

# remove the temporary file
rm "$temp_file"