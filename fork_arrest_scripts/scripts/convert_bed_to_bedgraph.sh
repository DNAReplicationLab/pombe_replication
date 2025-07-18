#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J convertBedToBedgraph
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Convert a bed file to bedgraph(s) using score/coverage for the fourth column

# usage
#------
# bash convert_bed_to_bedgraph.sh [-c] <input_bed_file> <fasta_fai_file> <output_prefix>
# can use sbatch if desired
# -c (optional): if used, discard bed file score column (if present) and use coverage of bed file instead
# input_bed_file: path to bed file
# fasta_fai_file: path to fasta fai file
# output_prefix: prefix for output bedgraph file(s)

# outputs
# -------
# bedgraph file(s) with the following naming convention:
# If a 3 column bed file is supplied and irrespective of whether -c is used or not,
#     <output_prefix>.coverage.bedgraph whose fourth column is coverage.
#     This is because 3 column bed files have no score data.
# If a bed file with at least six columns is supplied,
#     if -c is used, <output_prefix>.coverage.<strand>.bedgraph and <output_prefix>.coverage.bedgraph
#     if -c is not used, <output_prefix>.sum.<strand>.bedgraph and <output_prefix>.sum.bedgraph
#    where <strand> is either plus or minus and the additional file is all strands combined.
# If bed files have overlapping intervals, the scores are summed.

# Introduction to logic
# ---------------------
# There are many ways of converting a bed file to a bedgraph.
# The primary challenge is that a bedgraph generally contains one unique signal value per base in the genome,
# (and has no strand information) whereas bed files can contain multiple intervals that overlap the same base.
# Q: So, how do we collapse multiple intervals into a single value?
# A: See below.
# - If the input is a BED3 file, then there is no score or strand column, so we can only calculate coverage
#   i.e. how many times a base is overlapped by a bed interval irrespective of strand.
# - If the input is a BED6 file, then we can calculate coverage or we can aggregate the scores and we produce
#   three bedgraph files, one for each strand and one for all strands combined.
# - If the input is of some other format like BED3+1, then convert it to BED3 or BED6 before using this script.

# Details of logic
# ----------------
# For coverage calculation, the logic is straightforward.
# * At every base in the genome, calculate how many bed intervals overlap that base in a strand-specific manner
#   or otherwise depending on the input file.
# * Then, report this per-base measurement in a bedgraph file, combining adjacent bases with the same coverage
#   into a single interval.
# If we are aggregating scores, then the logic is similar, but we have made a choice to sum the scores.
# * At every base in the genome, we identify all the bed intervals that overlap that base
#   in a strand-specific manner or otherwise depending on the input file, then we can do one of the following:
#   (1) we can sum the scores
#   (2) we can calculate a per-base score per bed interval and then sum these per-base scores across all intervals
#   (3) we can do (1) or (2) but with some other function e.g. mean, median etc.
#   We have chosen to do number (1).
# * Then, report this per-base measurement in a bedgraph file, combining adjacent bases with the same value
#   into a single interval.

# Example input (treat spaces as tabs)
# -----------------------------------
# chrI 10 15 b 4 +
# chrI 15 20 b 5 +
# chrI 35 45 b 10 +
# chrI 40 50 b 11 +
# chrI 40 45 b 100 +
# chrI 40 46 b 200 +
# chrII 10 15 b 6 +

# Example output in the plus.bedgraph file if -c was not used
# -----------------------------------------------------------
# chrI 10 15 4
# chrI 15 20 5
# chrI 35 40 10
# chrI 40 45 321
# chrI 45 46 211
# chrI 46 50 11
# chrII 10 15 6

# Example output in the plus.bedgraph file if -c was used (there are more rows with 0 in them which are not shown)
# ----------------------------------------------------------------------------------------------------------------
# chrI 10 20 1
# chrI 35 40 1
# chrI 40 45 4
# chrI 45 46 2
# chrI 46 50 1
# chrII 10 15 1

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -python

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

# parse options
# -------------
# set defaults
discard_score=false

# parse arguments
while getopts ":c" opt; do
  case $opt in
    c)
      discard_score=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# parse arguments
# ---------------
# shift arguments away so that "$@" only contains positional arguments
shift "$((OPTIND-1))"

# assign arguments to variables
input_bed_file=${1:-}
fasta_fai_file=${2:-}
output_prefix=${3:-}

# check that the correct number of arguments were provided
if [ "$#" -ne 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 [-c] <input_bed_file> <fasta_fai_file> <output_prefix>"
    >&2 echo "input_bed_file: path to bed file"
    >&2 echo "fasta_fai_file: path to fasta fai file"
    >&2 echo "output_prefix: prefix for output bedgraph file(s)"
    >&2 echo "[-c] means the parameter is optional and if used, discard bed file score column (if present) and use coverage of bed file instead"
    exit 1;
fi

# check that the files exist
if [ ! -f "$input_bed_file" ]; then
    >&2 echo "ERROR: input bed file does not exist"
    exit 1;
fi

if [ ! -f "$fasta_fai_file" ]; then
    >&2 echo "ERROR: fasta fai file does not exist"
    exit 1;
fi

# make an output directory if needed
output_dir=$(dirname "$output_prefix")
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# sort the bed file
# -----------------
tmp_sorted_bed="$tmpDir"/sorted.bed
grep -E -v '^#|^browser|^track' "$input_bed_file" | sort -k 1,1 -k2,2n  > "$tmp_sorted_bed" 

# check if bed file is valid and then convert to bedgraph depending on whether score is to be discarded or not
# ------------------------------------------------------------------------------------------------------------
if [ ! "$(< "$tmp_sorted_bed" python validate_bed_against_fai.py "$fasta_fai_file" )" == "valid"  ]; then
  >&2 echo "Error: $input_bed_file does not have valid coordinates."
  exit 1;
fi

if [ "$(< "$tmp_sorted_bed" python validate_bed_format.py --six-columns --allow-float-score)" == "valid" ];
then
  if [ "$discard_score" = true ]; then
    bedtools genomecov -bga -g "$fasta_fai_file" -i "$tmp_sorted_bed" -strand + |\
      awk -v OFS=' ' '{print $1,$2,$3,$4}' > "$output_prefix".coverage.plus.bedgraph
    bedtools genomecov -bga -g "$fasta_fai_file" -i "$tmp_sorted_bed" -strand - |\
      awk -v OFS=' ' '{print $1,$2,$3,$4}' > "$output_prefix".coverage.minus.bedgraph
    bedtools genomecov -bga -g "$fasta_fai_file" -i "$tmp_sorted_bed"  |\
      awk -v OFS=' ' '{print $1,$2,$3,$4}' > "$output_prefix".coverage.bedgraph
  else
    < "$tmp_sorted_bed" grep -E -v "^#|^track|^browser" | awk -v OFS=" " '{if($6=="+"){print $1,$2,$3,$5}}' |\
      bash merge_overlapping_bedgraph_intervals.sh "$fasta_fai_file" > "$output_prefix".plus.bedgraph
    < "$tmp_sorted_bed" grep -E -v "^#|^track|^browser" | awk -v OFS=" " '{if($6=="-"){print $1,$2,$3,$5}}' |\
      bash merge_overlapping_bedgraph_intervals.sh "$fasta_fai_file" > "$output_prefix".minus.bedgraph
    < "$tmp_sorted_bed" grep -E -v "^#|^track|^browser" | awk -v OFS=" " '{print $1,$2,$3,$5}' |\
      bash merge_overlapping_bedgraph_intervals.sh "$fasta_fai_file" > "$output_prefix".bedgraph
  fi
elif [ "$(< "$tmp_sorted_bed" python validate_bed_format.py)" == "valid" ];
then
  bedtools genomecov -bga -g "$fasta_fai_file" -i "$tmp_sorted_bed" |\
    awk -v OFS=' ' '{print $1,$2,$3,$4}' > "$output_prefix".coverage.bedgraph
else
  >&2 echo "ERROR: input bed file is not valid"
  exit 1;
fi

# remove temporary directory
rm -rf "$tmpDir"