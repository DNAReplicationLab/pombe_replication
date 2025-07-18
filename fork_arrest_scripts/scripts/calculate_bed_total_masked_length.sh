#!/usr/bin/env bash

# goal
# =====
# Given a bed file, and a "mask" bed file, calculate the total length of the mask file
# that overlaps with the bed file.

# usage
# =====
# bash calculate_bed_total_masked_length.sh [-s|-S] <bed file> <mask bed file>
# NOTE: <> indicates required arguments which must be provided in the order shown.
#       Please fill in the path to the files you want to use in place of the placeholders (e.g. <bed file>).
# NOTE: [] indicates optional arguments.
#       -s: if used, intersect only if same strand.
#       -S: if used, intersect only if opposite strand.
#       if none of these options are used, then the strand information is ignored.
# NOTE: to use strand information, the bed file must have a 6th column with strand information.
# bed file: bed file of interest
# mask bed file: bed file of regions to mask

# example
# =======
# bash calculate_bed_total_masked_length.sh /path/to/bed/file /path/to/mask/bed/file
# bash calculate_bed_total_masked_length.sh -s /path/to/bed/file /path/to/mask/bed/file
# bash calculate_bed_total_masked_length.sh -S /path/to/bed/file /path/to/mask/bed/file

# output
# ======
# The first line is the total length of the mask bed file that overlaps with the bed file.

# fail on error
set -Eeuo pipefail

# process optional arguments
small_s=false
capital_s=false

while getopts "sS" opt; do
    case $opt in
        s)
            small_s=true
            ;;
        S)
            capital_s=true
            ;;
        \?)
            >&2 echo "ERROR: invalid option provided"
            exit 1
            ;;
    esac
done

# if both -s and -S are set, exit
if [ "$small_s" == "true" ] && [ "$capital_s" == "true" ]; then
    >&2 echo "ERROR: both -s and -S options are set. Please use only one of them."
    exit 1
fi

# shift arguments to exclude parsed options
shift $((OPTIND-1))

# check that the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash calculate_bed_total_masked_length.sh <bed file> <mask bed file>"
    exit 1
fi

# get input arguments
bed_file=$1
mask_bed_file=$2

# check that these files exist
if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file does not exist"
    exit 1
fi

if [ ! -f "$mask_bed_file" ]; then
    >&2 echo "ERROR: mask bed file does not exist"
    exit 1
fi

# load bedtools
source load_package.sh -bedtools -python

# depending on strand specific usages, check that the bed file has the right columns
if [ "$small_s" == "true" ] || [ "$capital_s" == "true" ]; then

  if [ ! "$(< "$bed_file" python validate_bed_format.py --no-dot-strand --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least six valid columns with sixth column = +/- when using -s/-S."
    exit 1;
  fi

  if [ ! "$(< "$mask_bed_file" python validate_bed_format.py --no-dot-strand --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least six valid columns with sixth column = +/- when using -s/-S."
    exit 1;
  fi

else

  if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least three valid columns."
    exit 1;
  fi

  if [ ! "$(< "$mask_bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least three valid columns."
    exit 1;
  fi

fi

# perform calculation
{
  if [ "$small_s" == "true" ]; then
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -s -i stdin -c 6 -o distinct | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "blank", 1000, $4}' |\
        bedtools intersect -s -a stdin -b "$mask_bed_file" -wo
  elif [ "$capital_s" == "true" ]; then
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -s -i stdin -c 6 -o distinct | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "blank", 1000, $4}' |\
        bedtools intersect -S -a stdin -b "$mask_bed_file" -wo
  else
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -i stdin |\
        bedtools intersect -a stdin -b "$mask_bed_file" -wo
  fi
} | awk 'BEGIN{a=0}{a+=$NF}END{print a}'