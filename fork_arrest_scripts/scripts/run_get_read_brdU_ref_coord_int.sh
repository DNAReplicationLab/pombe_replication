#!/bin/bash

#SBATCH --mem-per-cpu=7G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J aggModBamTrep
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# goal
# ====
# given a list of genomic windows, identify all reads that overlap with each window and output
# an analogue density and number of modified and unmodified bases per read id if any part of
# the read overlaps with the window.

# input and output
# =================
# Input: a modBAM file,
#        a file with a list of genomic windows, with one externally-measured value per window.
#        format of this file is tab-separated values with no header,
#        and columns are contig, start, end, irrelevant, value.
#        contig, start, end refer to coordinates on the reference.
#        irrelevant is a column whose values we won't be using.
#        value is some numeric value.
# Output: Output is sent to the specified file.
#         Output format is tab-separated with no headers.
#         Columns are contig, start, end, irrelevant, value, read_id, analogue_mean, num_bases.
#         The first five are just a repeat of the input.
#         The last three are explained below.
#         For each window, we calculate and output all (read id, analogue density, number of mod + unmod bases) tuples
#         where read id refers to any read id that overlaps with the genomic window (even 1 base of overlap is
#         considered valid here), analogue density is the mean analogue density over the bases of the read in the
#         window, and number of mod + unmod bases refers to number of bases in the read in the genomic window.

# stop execution if any step fails
set -e

# get config directory info
source config.sh

# read command line inputs and set variables
if [ "$#" -ne 4 ]; then
    echo "Usage: bash <program_name.sh> mod_bam bed_file output_file temp_directory"
    exit;
fi

# set input and output directories, files
modBAM=$1
bedFile=$2
outFile=$3

# check that input files exist
if [ ! -f "$modBAM" ]; then
    echo "Error: modBAM file $modBAM does not exist."
    exit;
fi

if [ ! -f "$modBAM".bai ]; then
    echo "Error: modBAM file index $modBAM.bai does not exist."
    exit;
fi

if [ ! -f "$bedFile" ]; then
    echo "Error: bed file $bedFile does not exist."
    exit;
fi

# make temporary directory
mkdir -p "$4";
# shellcheck disable=SC2164
tmpDir=$(cd "$4"; pwd)

# load git repo labels
source load_git_repo_labels.sh
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# load python
source load_package.sh -python

# split bed file into smaller files w N lines per file (except last file)
# generate a random prefix for split file name
rndStr=${tmpDir}/$(openssl rand -hex 6)
split --verbose -l500 "$bedFile" "${rndStr}".

for piece in "${rndStr}".*;
do
    # calculate aggregate BrdU statistics per read per window
    < "$piece" \
        sed '1icontig start end dummy trep' |\
        python get_read_brdU_ref_coord_int.py \
            --modBAM "$modBAM" --thres 0.5 \
        > "${piece}".piece &
done

wait;

# concatenate results
cat "${rndStr}".*.piece |\
    sed "1i# using modBAM ${modBAM}" |\
    sed "1i# using Trep file ${bedFile}" |\
    sed "1i${infoStr}" \
    > "$outFile"

# delete temporary files
rm "${rndStr}".*;
