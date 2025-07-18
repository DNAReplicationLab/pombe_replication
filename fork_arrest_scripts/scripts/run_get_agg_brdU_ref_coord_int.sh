#!/bin/bash

#SBATCH --mem=200G
#SBATCH -c 1 
#SBATCH -p ei-medium
#SBATCH -J aggModBamTrep
#SBATCH --mail-type=END,FAIL
#SBATCH --time=47:59:59

# stop execution if any step fails
set -e

# introduction
# ==============
# calculate mean BrdU amounts and total number of modified and unmodified bases per genomic window.

# typical usage
# =============
# bash <program_name.sh> $modbam $bed $output
# can use sbatch in place of bash
# modbam: .mod.bam file with analogue modification probabilities
# bed: .bed file with five tab-separated columns and no header: contig, start, end, irrelevant, quantity
#       contig, start, end are coordinates on a reference genome.
#       irrelevant is an irrelevant column that we won't be using.
#       quantity is some quantity that arises out of a measurement and is associated with that genomic window.
# output: output file. has same format as input bed file but with two additional columns that report mean BrdU and
#         total number of modified and unmodified bases per window.

# check number of arguments
if [ "$#" -ne 3 ]
then
  echo "Incorrect number of arguments! "
  echo "Usage: sbatch <program_name.sh> modbam_file bed_file output_file"
  exit
fi

# set input and output directories, files
modBAM=$1
bedFile=$2
outFile=$3

# check that modBAM file exists
if [ ! -f "$modBAM" ]
then
  echo "modBAM file does not exist!"
  exit
fi

# check that bed file exists
if [ ! -f "$bedFile" ]
then
  echo "bed file does not exist!"
  exit
fi

# load git repo labels
source load_git_repo_labels.sh
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# load python
source load_package.sh -python

# calculate aggregate BrdU statistics per window
< "$bedFile" \
    sed '1icontig start end dummy trep' |\
    python get_agg_brdU_ref_coord_int.py \
        --modBAM "$modBAM"\
        --thres 0.5 |\
    sed "1i# using modBAM ${modBAM}" |\
    sed "1i# using Trep file ${bedFile}" |\
    sed "1i${infoStr}" \
        > "$outFile"