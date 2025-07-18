#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J convertBamToBedgraph
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Convert a bam file to bedgraph(s) using coverage of mapped reads for the fourth column

# usage
#------
# bash convert_bam_to_bedgraph.sh <input_bam_file> <fasta_fai_file> <output_prefix> [region]
# can use sbatch if desired
# input_bam_file: path to bam file
# fasta_fai_file: path to fasta fai file
# output_prefix: prefix for output bedgraph file(s)
# region (optional): region to restrict the analysis to reads mapping to that region. Format is "chr:start-end" or
#                    "chr". Default is the whole genome.

# outputs
# -------
# bedgraph file(s) with the following naming convention:
#    <output_prefix>.coverage.<strand>.bedgraph and <output_prefix>.coverage.bedgraph
#    where <strand> is either plus or minus and the additional file is all strands combined.

# Details of logic
# ----------------
# Convert bam to a bed file using only mapped reads, and then use convert_bed_to_bedgraph.sh.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# assign arguments to variables
input_bam_file=${1:-}
fasta_fai_file=${2:-}
output_prefix=${3:-}
region=${4:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 <input_bam_file> <fasta_fai_file> <output_prefix> [region]"
    >&2 echo "input_bam_file: path to bam file"
    >&2 echo "fasta_fai_file: path to fasta fai file"
    >&2 echo "output_prefix: prefix for output bedgraph file(s)"
    >&2 echo "region (optional): region to restrict the analysis to reads mapping to that region. "
    >&2 echo "                   Format is 'chr:start-end' or 'chr'. Default is the whole genome."
    exit 1;
fi

# check that the files exist
if [ ! -f "$input_bam_file" ] || [ ! -f "$input_bam_file.bai" ]; then
    >&2 echo "ERROR: input bam file or index does not exist"
    exit 1;
fi

if [ ! -f "$fasta_fai_file" ]; then
    >&2 echo "ERROR: input fasta fai file does not exist"
    exit 1;
fi

# check that region is of the correct format
if [ -n "$region" ]; then
    if [[ ! "$region" =~ ^[A-Za-z0-9_.]+:[0-9]+-[0-9]+$ ]] && [[ ! "$region" =~ ^[A-Za-z0-9_.]+$ ]]; then
        >&2 echo "ERROR: region is not of the correct format: $region"
        >&2 echo "Correct format is 'contig:start-end'"
        >&2 echo "NOTE: if your contig contains unusual characters, the script may not work"
        exit 1;
    fi
fi

# convert bam file to bed
tmp_bed_file="$tmpDir"/"$(basename "$input_bam_file" .bam)".bed
# shellcheck disable=SC2086
samtools view -h -b "$input_bam_file" --exclude-flags UNMAP $region |\
  bedtools bamtobed -i - > "$tmp_bed_file"
bash convert_bed_to_bedgraph.sh -c "$tmp_bed_file" "$fasta_fai_file" "$output_prefix"

# remove temporary directory
rm -rf "$tmpDir"