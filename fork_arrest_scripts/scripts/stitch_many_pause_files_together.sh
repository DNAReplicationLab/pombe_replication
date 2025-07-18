#!/bin/bash

#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J stitchPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --time 00:30:00

# goal
# -----
# Stitch together many pause files into one.
# Each pause file is in a CSV/CSV-like format with headers.

# usage
# ------
# sbatch stitch_many_pause_files_together.sh inputFilePrefix opFile
# inputFilePrefix: prefix of input files. Files with names inputFilePrefix_part_* are stitched together and deleted.
#                  * is a wildcard that matches any string of any length from 0 to infinity.
# opFile: output file

# stop execution if any command fails
set -e

# get inputs
input_file=$1
output_file=$2

# check that two input arguments have been supplied
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "usage: sbatch stitch_many_pause_files_together.sh input_file_prefix output_file"
    exit 1
fi

# load package
source load_package.sh -miller

# stitch together the files and delete the input files
mlr --csv --ifs ' ' --ofs ' ' --skip-comments cat "$input_file"_part_* > "$output_file";
rm -rf "$input_file"_part_*;