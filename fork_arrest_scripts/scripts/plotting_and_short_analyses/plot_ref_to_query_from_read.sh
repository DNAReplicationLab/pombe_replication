#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotRefToQueryFromRead
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1:59:59
#SBATCH --constraint=""

# goal
# -----
# Given a read id, plot XY scatter X = coordinate on reference, Y = coordinate on query

# usage
#------
# bash plot_ref_to_query_from_read.sh read_id bam_file output_prefix
# read_id: read id
# bam_file: bam file with alignments, must be indexed. NOTE: If you use the mod bam file obtained from a detect file,
#           you will get a perfect Y = X line as for convenience we store modification data assuming a perfect match
#           between read and reference as all our modification calls are on the reference.
# output_prefix: prefix for output files. e.g. If you want files with names starting with blah under the directory
#                /path/to/dir, then set output_prefix to /path/to/dir/blah

# warning
# -------
# if there are multiple reads with the same read_id, the script will fail.

# outputs
# -------
# A tab-separated text file and plot in png format (output_prefix_{read_id}.tsv and output_prefix_{read_id}.png)
# where {read_id} is replaced with the read id provided as input.
# For what's in the text file, consult the called python script below.

# stop execution if any command fails
set -e

# set directory paths
curr_dir=$(pwd)
main_dir=$(cd ..; pwd)

# change to scripts directory
cd "$main_dir"

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -samtools -python -miller -R

# load configuration
source config.sh

# change back to original directory
cd "$curr_dir"

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# assign arguments to variables
read_id=$1
bam_file=$2
output_prefix=$(realpath "$3")

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 read_id bam_file output_prefix"
    >&2 echo "read_id: read id"
    >&2 echo "bam_file: bam file with alignments, must be indexed. NOTE: do not use mod bam file obtained from a detect file here"
    >&2 echo "output_prefix: prefix for output files"
    >&2 echo "For more details, consult the comments at the top of the script"
    >&2 echo "Exiting."
    exit 1;
fi

# check that bam file exists
if [ ! -f "$bam_file" ]; then
    >&2 echo "ERROR: bam file provided does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

# check that bam file is indexed
if [ ! -f "$bam_file.bai" ]; then
    >&2 echo "ERROR: bam file provided is not indexed"
    >&2 echo "Exiting."
    exit 1;
fi

# extract read of interest
samtools view -e 'qname=="'"$read_id"'"' "$bam_file" > "$tmpDir/$read_id.sam"

# check that only one read was found, otherwise exit with message
if [ "$(wc -l < "$tmpDir/$read_id.sam")" -ne 1 ]; then
    >&2 echo "ERROR: read id provided does not uniquely identify a read in the bam file provided"
    >&2 echo "Exiting."
    exit 1;
fi

# extract cigar string and output text file with X, Y information
cd "$main_dir"
< "$tmpDir/$read_id.sam"  awk -v OFS="\t" '{print $1,$4,$6}' |\
  sed '1iindex\tposOnRefStart\tcigarStr' | python get_ref_to_query_from_cigar.py |\
  awk -F" " -v OFS="\t" '{print $1,$2/1000,$3/1000}' | sed '1iread_id\tposOnRef\tposOnRead' \
  > "${output_prefix}_${read_id}.tsv"
cd "$curr_dir"

# extract min and max of X and Y
x_min=$(mlr --tsv --headerless-csv-output stats1 -f posOnRef -a min "${output_prefix}_${read_id}.tsv")
x_max=$(mlr --tsv --headerless-csv-output stats1 -f posOnRef -a max "${output_prefix}_${read_id}.tsv")
y_min=$(mlr --tsv --headerless-csv-output stats1 -f posOnRead -a min "${output_prefix}_${read_id}.tsv")
y_max=$(mlr --tsv --headerless-csv-output stats1 -f posOnRead -a max "${output_prefix}_${read_id}.tsv")

# make the plot
< "${output_prefix}_${read_id}.tsv" Rscript plot_scatter.R posOnRef posOnRead "${output_prefix}_${read_id}.png"\
  "Reference coordinate (kb)" "Read coordinate (kb)" "$x_min,$x_max" "$y_min,$y_max" 1 1 0 0 18 12 50