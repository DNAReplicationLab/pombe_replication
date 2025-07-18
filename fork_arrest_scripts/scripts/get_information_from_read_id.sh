#!/bin/bash

# Goal
# ====
# This script takes a bam file and a read id and returns the alignment information for that read id

# usage
# =====
# bash get_information_from_read_id.sh [-o <bam_file_with_just_read>] <bam file> <read id>
# bam file: bam file to be used
# read id: read id to be used
# -o: optional argument to create the specified output bam file with just the read id
# To specify this, you must remove the [] and write the stuff inside the [].

# output
# ======
# output is to stdout and looks like (some of the spaces are tabs, the lines start with #, and order of the lines
# may change):
## bam: <bam file>
## read id: <read id>
## alignment information:  <chromosome> <start> <end> <read id> <score> <+/->
## orientation: <+/->

# fail if any command fails
set -e

# process flags
output_bam_file=/dev/null
while getopts ":o:" opt; do
  case $opt in
    o) output_bam_file="$OPTARG"
    ;;
    \?) : # do nothing
    ;;
  esac
done

# shift arguments to exclude parsed options
shift $((OPTIND-1))

print_usage() {
  >&2 echo "$1"
  >&2 echo "usage: bash get_information_from_read_id.sh [-o <bam_file_with_just_read>] <bam file> <read id>"
  >&2 echo "NOTE: [] means the argument is optional. If you want to specify this, then you must remove the []"
  >&2 echo "The optional argument -o can be used to make an output bam file with just the read id"
  exit 1
}

# check if the correct number of arguments is given
if [ "$#" -ne 2 ]; then
  print_usage "Illegal number of parameters"
fi

# check if the bam file exists
if [ ! -f "$1" ]; then
  print_usage "bam file does not exist"
fi

# load bedtools, samtools
source load_package.sh -bedtools -samtools -python -modkit;

# make a temporary file to store bam information
temp_bam_file=$(mktemp)

# print output
# shellcheck disable=SC2015
samtools view -b -e 'qname=~"'"$2"'"' "$1" > "$temp_bam_file"
count_entries=$(samtools view -c "$temp_bam_file")

# if there are no alignments, delete temporary file and exit with error
# shellcheck disable=SC2046
if [ "$count_entries" -eq 0 ]; then
  >&2 echo "# no alignment information found";
  rm "$temp_bam_file";
  exit 1;
fi

# if there are multiple alignments, print a warning
# shellcheck disable=SC2046
if [ "$count_entries" -gt 1 ]; then
  echo "# warning: multiple alignments found";
fi

# print alignment information
bedtools bamtobed -i "$temp_bam_file" |\
  awk '{print "# alignment information: ", $0, "\n# orientation:", $6}'

# print bam file and read id
echo "# bam:" "$1"
echo "# read id:" "$2"

# get custom tags
samtools view -h "$temp_bam_file" | python reformat_custom_bam_annotation_tags.py

# print mod summary
echo "# Modkit summary follows."
echo "# Please note this is for a quick summary; check modkit documentation to see if ref-dependent "
echo "# or ref-independent coords are used and what thresholds are used if you want to use these numbers further."
modkit summary --no-sampling --no-filtering "$temp_bam_file" | sed 's/^/#/'

# if output bam file is specified, sort the temporary bam file and save it
if [ "$output_bam_file" != "/dev/null" ]; then
  samtools sort -o "$output_bam_file" "$temp_bam_file"
  samtools index "$output_bam_file"
fi

# delete temporary file
rm "$temp_bam_file"