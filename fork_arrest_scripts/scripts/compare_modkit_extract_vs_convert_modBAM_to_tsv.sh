#!/bin/bash

#SBATCH --mem-per-cpu=40G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J compareModOutput
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Compare the output of modkit extract and convert_modBAM_to_tsv.py.
# We want not even a single line to be different between the two outputs.
# Due to historical reasons, our codebase requires a python script to read
# modBAM files. ONT used to support this through a package called modbampy
# but have unfortunately discontinued it and replaced it with a command line
# tool called modkit which is somewhat inconvenient for us to use as it
# is not a python package. Thus, we have written our own modBAM parser
# and are testing it by comparing its output with modkit's output.

# usage
#------
# sbatch compare_modkit_extract_vs_convert_modBAM_to_tsv.sh $mod_bam_file $base $mod [sample_frac]
# NOTE: can also use bash if you want to run it locally
# NOTE: [] denotes optional arguments. To use optional arguments, all previous optional arguments must be provided,
#       and the [] must be removed.
# mod_bam_file: the modBAM file with modification information. NOTE: the file must not have multiple alignments per
#               the same readID.
# base: the base to extract
# mod: the modification code to extract
# sample_frac: (default 1) the fraction of reads to sample. NOTE that if the modBAM file is large, you will
#              run into memory issues if you do not sample. We recommend a read number in the hundreds.
#              If you want to set this to 1, please set this to a number very close to 1 but not 1,
#              e.g. 0.999999.

# outputs
# -------
# A few lines are printed to stdout showing the number of lines that do not exist in both tables,
# and the maximum of the square of the difference in modification quality and reference position per base
# between the two tables.
# Ideally, the first number should be 0, and the second pair of numbers should be 0 (we cannot expect
# modification probability differences to be exactly 0 due to floating point errors and due to
# the 1/256 quantization of the modification quality scores in the modBAM file).

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -samtools -modkit -miller

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# assign arguments to variables
input_mod_bam=${1:-}
base=${2:-}
mod=${3:-}
sample_frac=${4:-0.9999999999999}

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch compare_modkit_extract_vs_convert_modBAM_to_tsv.sh \$mod_bam_file \$base \$mod [sample_frac]"
    >&2 echo "mod_bam_file: the modBAM file with modification information, must not have multiple alignments per read id"
    >&2 echo "base: the base to extract"
    >&2 echo "mod: the modification code to extract"
    >&2 echo "sample_frac: (default 1) the fraction of reads to sample. To avoid memory issues, we recommend setting "
    >&2 echo "             this such that the resultant number of reads is in the hundreds. To set to 1, please use "
    >&2 echo "             a number very close to 1 but not 1 e.g. 0.9999 "
    exit 1;
fi

# set up some output tsv temporary files
output_tsv_folder=$tmpDir/output_tsv
mkdir -p "$output_tsv_folder"
output_tsv_1=$output_tsv_folder/modkit_extract.tsv
output_tsv_2=$output_tsv_folder/convert_modBAM_to_tsv.tsv
output_tsv_3=$output_tsv_folder/merged.tsv
output_tsv_4=$output_tsv_folder/diff.tsv

# subsample the modBAM file using samtools
samtools view -s "$sample_frac" -b "$input_mod_bam" > "$tmpDir/subsampled_modBAM.bam"
samtools sort -o "$tmpDir/subsampled_modBAM_sorted.bam" "$tmpDir/subsampled_modBAM.bam"
samtools index "$tmpDir/subsampled_modBAM_sorted.bam"

# first run modkit extract
modkit extract full "$tmpDir/subsampled_modBAM_sorted.bam" "$output_tsv_1".tmp
# shellcheck disable=SC2016
mlr --tsv filter '$canonical_base =="'"$base"'" && $modified_primary_base == "'"$mod"'"' \
  "$output_tsv_1".tmp > "$output_tsv_1"

# then use our own script to perform the conversion
samtools view -h "$tmpDir/subsampled_modBAM_sorted.bam" |\
  python convert_modBAM_to_tsv.py --tag "$mod" --base "$base" > "$output_tsv_2"

# merge the two files
mlr --tsv --skip-comments join -f "$output_tsv_1" -j read_id,forward_read_position --lp "mk_" --rp "us_" \
  "$output_tsv_2" > "$output_tsv_3"

# get the difference
mlr --tsv --skip-comments join -f "$output_tsv_1" -j read_id,forward_read_position --lp "mk_" --rp "us_" \
  --np --ul --ur "$output_tsv_2" > "$output_tsv_4"

# check output_tsv_4 has no lines
echo "Number of lines that do not exist in both tables (we'd prefer this to be zero): $(wc -l "$output_tsv_4")"

echo "Maximum of the square of the difference in modification quality and "
echo "in position along the reference genome per base between the two tables: "
# find max of difference
# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --tsv put '$sq_diff=($mk_mod_qual - $us_mod_qual)*($mk_mod_qual - $us_mod_qual);'\
'$ref_diff=($mk_ref_position - $us_ref_position)*($mk_ref_position - $us_ref_position);' \
  then stats1 -a max -f sq_diff,ref_diff "$output_tsv_3"

echo "Please note that we want both these quantities to be zero, but the modification "
echo "quality difference may not be identically zero due to the 1/256 quantization "
echo "required by the modBAM format."

# remove temporary files
rm -rf "$tmpDir"
