#!/bin/bash
#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J joinBasecalledLenAlignLen
#SBATCH --mail-type=END,FAIL
#SBATCH --time=01:00:00

# Goal
# ----
# Every read id has a corresponding alignment length(s) and basecalled length, so we tabulate them here.

# Usage
# -----
# sbatch join_alignLen_basecalledLen_per_readid.sh sequencing_summary bam_file.bam output_file
# sequencing_summary: sequencing summary output by the basecaller guppy.
# bam_file.bam: Any bam file (output by minimap2, a mod bam file etc.) with some reads from the sequencing_summary file.
# output_file: Output file name.

# Output
# ------
# output_file will contain three tab-separated columns of data with header: read_id, align_length,
# sequence_length_template. All lengths are in bases, not kb.
# NOTE: if bam file contains multiple alignments for a read id, then there will be multiple lines for that read id
# in the output file.
# NOTE: we do not anticipate the sequencing_summary file to contain multiple lines for a read id.

# File format note
# ----------------
# sequencing_summary is a tab delimited file with header. we'll be using the columns read_id and
# sequence_length_template which contains the basecalled length in b.
# bam file just follows the standard bam format.

# fail on error
set -e

# load bedtools, miller
source load_package.sh -bedtools -miller

# check that there are three arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: sbatch join_alignLen_basecalledLen_per_readid.sh sequencing_summary bam_file.bam output_file"
  exit 1
fi

# check that the sequencing summary file exists
if [ ! -f "$1" ]; then
  echo "Error: sequencing summary file $1 does not exist"
  exit 1
fi

# check that the bam file exists
if [ ! -f "$2" ]; then
  echo "Error: bam file $2 does not exist"
  exit 1
fi

# check that the output file does not exist
if [ -f "$3" ]; then
  echo "Error: output file $3 already exists"
  exit 1
fi

# convert bam to bed and send to a temporary file
tmp_file=$(mktemp)
bedtools bamtobed -i "$2"  | awk '{print $4 "\t" $3 - $2}' | sed '1iread_id\talign_length' > "$tmp_file"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
{
  echo "# input file sequencing summary: $1"
  echo "# input file bam: $2"
  mlr --itsv --otsv join -j read_id -f "$tmp_file" \
      then cut -f read_id,align_length,sequence_length_template "$1";
} > "$3"

# remove temporary file
rm "$tmp_file"