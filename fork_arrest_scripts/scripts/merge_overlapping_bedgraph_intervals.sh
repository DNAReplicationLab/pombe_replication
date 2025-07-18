#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J mergeOverlappingBedgraphIntervals
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Merge overlapping intervals in a bedgraph file, summing the values. Non-overlapping intervals are not affected.
# Bedgraphs are supposed to contain non-overlapping intervals as they are a record of a signal vs genome coordinate
# (where coordinates are expressed in windows). When intervals overlap, the signal is ambiguous.

# usage
#------
# cat <input_bedgraph_file> | bash merge_overlapping_bedgraph_intervals.sh $fasta_fai > <output_bedgraph_file>
# can use sbatch if desired
# fasta_fai: path to fasta fai file

# outputs
# -------
# Output bedgraph where overlapping intervals have been merged and values summed.

# Example input
# -------------
# chrI 10 15 4
# chrI 15 20 5
# chrI 35 45 10
# chrI 40 50 11
# chrI 40 45 100
# chrI 40 46 200
# chrII 10 15 6

# Example output
# --------------
# chrI 10 15 4
# chrI 15 20 5
# chrI 35 40 10
# chrI 40 45 321
# chrI 45 46 211
# chrI 46 50 11
# chrII 10 15 6

# Example logic
# -------------
# chrI 40 - 45 was a common region. So, we summed the values 10, 11, 100, 200 to get 321.
# chrI 45 - 46 was a common region. So, we summed the values 11, 200 to get 211.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -python -bedops

# load configuration
source config.sh

# set temporary directory
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

# check that the correct number of arguments were provided
if [ "$#" -ne 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: cat <input_bedgraph_file> | bash merge_overlapping_bedgraph_intervals.sh <fasta_fai_file> > <output_bedgraph_file>"
    >&2 echo "fasta_fai_file: path to fasta fai file"
    exit 1;
fi

fasta_fai_file=${1:-}

# check that input file exists
if [ ! -f "$fasta_fai_file" ]; then
    >&2 echo "ERROR: fasta fai file does not exist"
    exit 1;
fi

# create temporary files
input_bed="$tmpDir"/input.bed
coverage_bed="$tmpDir"/coverage.bed

# read input bedgraph file from stdin and send to temporary file in bed3+1 format
grep -E -v '^track|^browser|^#' /dev/stdin | awk -v OFS='\t' '{if($2!=$3){print $1,$2,$3,"blank",$4}}' |\
  sort -k 1,1 -k2,2n  > "$input_bed"

# if input bedgraph file is empty, exit
if [ "$(wc -l < "$input_bed")" -eq 0 ]; then
  rm -r "$tmpDir"
  exit 0;
fi

# check that the bed file is valid
if [ ! "$(< "$input_bed" python validate_bed_against_fai.py "$fasta_fai_file" )" == "valid"  ]; then
  >&2 echo "Error: input bed does not have valid coordinates."
  exit 1;
fi

if [ ! "$(< "$input_bed" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: input bed does not have valid data."
  exit 1;
fi

# This step identifies windows that overlap and produces windows with adjusted boundaries.
bedops --partition "$input_bed" | sort -k 1,1 -k2,2n  > "$coverage_bed"

# Sum values over windows identified in previous step
bedtools map -a "$coverage_bed" -b "$input_bed" -c 5 -o sum |\
 awk -v OFS=' ' '{print $1,$2,$3,$4}' |\
 insert_calling_script_header "$@"

# remove temporary directory
rm -r "$tmpDir"