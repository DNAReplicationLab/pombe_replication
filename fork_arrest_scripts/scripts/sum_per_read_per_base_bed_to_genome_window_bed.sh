#!/bin/bash

#SBATCH --mem-per-cpu=100G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J normPerReadPerBaseBedToGenomeBed
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# given a bed file with a per-read per-base signal, construct genomic windows, and sum intersecting values per window.

# usage
#------
# bash sum_per_read_per_base_bed_to_genome_window_bed.sh input.bed input.fai window_size
# input.bed: bed file with per-read per-base signal
# input.fai: fasta index file for the genome, generate with samtools faidx fasta_file if not available
# window_size: size of the window to tile the genome in bp

# outputs
# -------
# to standard output in bedgraph format: i.e. space-separated with no column names with
# the columns contig, start, end, sum
# where contig, start, end are window boundaries and sum is the sum of the values in the window.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -bedtools -python

# load configuration
source config.sh

# set temporary directory and make it
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
input_bed=$1
input_fai=$2
window_size=$3

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
  >&2 echo "ERROR: incorrect number of arguments provided"
  >&2 echo "usage: bash sum_per_read_per_base_bed_to_genome_window_bed.sh input.bed input.fai window_size"
  >&2 echo "input.bed: bed file with per-read per-base signal"
  >&2 echo "input.fai: fasta index file for the genome"
  >&2 echo "window_size: size of the window to tile the genome in bp"
  exit 1;
fi

# if bed_file is stdin or -, then read from stdin, store in a temporary file, and set input_bed to that temporary file
if [ "$input_bed" == "stdin" ] || [ "$input_bed" == "-" ]; then
  tmp_bed_file_stdin=$(mktemp -p "$tmpDir" tmp_bed_file_stdin_XXXXXX.bed)
  cat > "$tmp_bed_file_stdin"
  input_bed="$tmp_bed_file_stdin"
fi

# check that the input bed file exists
if [ ! -f "$input_bed" ]; then
  >&2 echo "ERROR: input bed file $input_bed does not exist"
  exit 1;
fi

# check that the input fai file exists
if [ ! -f "$input_fai" ]; then
  >&2 echo "ERROR: input fai file $input_fai does not exist"
  exit 1;
fi

# check that the window size is a positive integer
if ! [[ "$window_size" =~ ^[0-9]+$ ]]; then
  >&2 echo "ERROR: window size $window_size is not a positive integer"
  exit 1;
fi

if [ "$window_size" -lt 1 ]; then
  >&2 echo "ERROR: window size $window_size is not a positive integer"
  exit 1;
fi

# validate bed format
if [ ! "$(< "$input_bed" python validate_bed_format.py --six-columns --allow-float-score\
                --require-uuid --single-base --score-restrict-to-0-0pt05)" == "valid" ];
then
  >&2 echo "Error: bed file is invalid. please have a look at the script for the criteria we use to check this."
  exit 1;
fi

if [ ! "$(< "$input_bed" python validate_bed_against_fai.py "$input_fai")" == "valid" ];
then
  >&2 echo "Error: bed file does not match fai file."
  exit 1;
fi

# sort bed file
input_bed_sorted=$(mktemp -p "$tmpDir" tmp_input_bed_sorted_XXXXXX.bed)
bedtools sort -faidx "$input_fai" -i "$input_bed" > "$input_bed_sorted"

# perform the calculation
# NOTE: (1) we don't have to think about overlap between generated windows and windows in our input_bed_sorted data
#           in the bedtools map command as the bed file contains intervals of size 1 bp that are non-overlapping.
#       (2) makewindows makes a mis-sized window at the end of each contig if the contig length is not a
#           multiple of window_size; this is alright I think.
bedtools makewindows -g "$input_fai" -w "$window_size" |\
    bedtools map -a - -b "$input_bed_sorted" -c 5 -o sum -null 0 |\
    awk -v OFS=" " '{print $1,$2,$3,$4}' |\
    insert_calling_script_header "$@"

# remove temporary files
rm -rf "$tmpDir"