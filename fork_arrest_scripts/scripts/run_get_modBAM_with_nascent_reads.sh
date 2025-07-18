#!/bin/bash
#SBATCH --mem-per-cpu=7G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J sbsModBamNascent
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# goal
# ====
# given a mod bam file, extract the reads that are nascent into a new bam file and sort, index it.

# usage, input, and output
# ========================
# usage: sbatch run_get_modBAM_with_nascent_reads.sh <input> <output> <temporary_directory>
# input: a modBAM file
# output: output is sent to the specified bam file name.
# temporary_directory: a directory to store temporary files. Must exist already and
#                      will not be deleted at end of script.

# stop execution if any step fails
set -e

# get config directory info
source config.sh

# read command line inputs and set variables
if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch <program_name.sh> input_mod_bam_file output_mod_bam_file tmp_dir"
    echo "tmp_dir and input_mod_bam_file must exist."
    exit;
fi

# set input and output files
modBAMInput=$1
modBAMOutput=$2
tmpDir=$3

# check that input file exists
if [ ! -f "$modBAMInput" ]; then
    echo "Input file $modBAMInput does not exist."
    exit;
fi

# check that the temp directory exists
if [ ! -d "$tmpDir" ]; then
    echo "Temp directory $tmpDir does not exist."
    exit;
fi

# load python
source load_package.sh -python -samtools

# create a temporary file
tmpFile=$(mktemp -p "$tmpDir")

# extract the reads that are nascent to the temporary file
samtools view -h "$modBAMInput" |\
  python get_modBAM_with_nascent_reads.py --window 300 |\
  samtools sort - |\
  samtools view -Sb -o "$tmpFile" -

# sort output into the out file
samtools sort -@ 16 -o "$modBAMOutput" -T "$tmpDir" "$tmpFile"

# index the file
samtools index "$modBAMOutput"

# remove the temporary file
rm "$tmpFile"