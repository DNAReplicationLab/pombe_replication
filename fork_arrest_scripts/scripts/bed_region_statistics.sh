#!/bin/bash

# goal
# -----
# print some statistics about a bed file

# usage
#------
# bash bed_region_statistics.sh <bed file> <prefix> <fasta_file>
# bed file: bed file to print statistics about
# prefix: (optional) this word and an underscore are attached to each measurement name in the output json.
#         not used if not set.
# fasta_file: (optional) fasta file to calculate base content, which is not output if not set.

# outputs
# -------
# prints statistics about the bed file to stdout

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -python -miller -jq > /dev/null 2>&1

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# assign arguments to variables
bed_file=${1:-}
prefix=${2:-}
fasta_file=${3:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 <bed file> <prefix> <fasta_file>"
    >&2 echo "bed file: bed file to print statistics about"
    >&2 echo "prefix: (optional) this word and an underscore are attached to each measurement name in the output json."
    >&2 echo "        not used if not set."
    >&2 echo "fasta_file: (optional) fasta file to calculate base content, which is not output if not set."
    exit 1;
fi

# if bed_file is stdin, then read from stdin, store in a temporary file, and set bed_file to that temporary file
if [ "$bed_file" == "stdin" ]; then
  tmp_bed_file_stdin=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  cat > "$tmp_bed_file_stdin"
  bed_file="$tmp_bed_file_stdin"
fi

# check that the input file exists
if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: input file does not exist"
    exit 1;
fi

# check that the input bed file is valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --six-columns --allow-float-score --no-dot-strand)" == "valid" ];
then
  >&2 echo "Error: bed file needs at least six columns and must not have dot strand."
  exit 1;
fi

# if the fasta file exists, then validate the bed file against the fasta .fai file
if [ -f "$fasta_file" ]; then

  if [ ! -f "$fasta_file".fai ]; then
    >&2 echo "Error: fasta index file $fasta_file.fai does not exist."
    exit 1;
  fi

  # bed file must have contigs only in the fasta index file
  if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fasta_file".fai)" == "valid" ]; then
    >&2 echo "Error: bed file must have contigs only in the fasta index file."
    exit 1;
  fi

fi

# if prefix is not null, then attach an underscore to it
if [ -n "$prefix" ]; then
  prefix="${prefix}_"
fi

# extract first six columns from the bed file
tmp_bed_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
< "$bed_file" grep -E -v '^browser|^track|^#' |\
  sort -k 1,1 -k2,2n |\
  awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1,$2,$3,$4,$5,$6}}' > "$tmp_bed_file"

# merge file along each strand
tmp_merged_bed_file=$(mktemp -p "${config[scratchDir]:-}"/tmp);
bedtools merge -s -i "$tmp_bed_file" -c 6 -o distinct |\
  awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "blank", 1000, $4}' > "$tmp_merged_bed_file"

# merge bed file ignoring strand
tmp_merged_bed_file_ignore_strand=$(mktemp -p "${config[scratchDir]:-}"/tmp);
bedtools merge -i "$tmp_bed_file" |\
  awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > "$tmp_merged_bed_file_ignore_strand"

# print outputs
echo "{"
if [ "$bed_file" == "$tmp_bed_file_stdin" ]; then
  echo "\"${prefix}bed_file\": \"stdin\","
else
  echo "\"${prefix}bed_file\": \"$bed_file\","
fi
echo "\"${prefix}number_of_regions\": $(wc -l < "$tmp_bed_file"),"
echo "\"${prefix}number_of_+_regions\": $(< "$tmp_bed_file" awk 'BEGIN{a=0}$6 == "+"{a+=1}END{print a}'),"
echo "\"${prefix}number_of_-_regions\": $(< "$tmp_bed_file" awk 'BEGIN{a=0}$6 == "-"{a+=1}END{print a}'),"
echo "\"${prefix}number_of_._regions\": $(< "$tmp_bed_file" awk 'BEGIN{a=0}$6 == "."{a+=1}END{print a}'),"
echo "\"${prefix}total_length_bp\": $(< "$tmp_bed_file" awk 'BEGIN{a=0}{a+=$3-$2}END{print a}'),"

# obtain number of intersecting intervals in ROI (excluding self-intersections)
# shellcheck disable=SC1010
# shellcheck disable=SC2016
twice_number_of_intersections=\
$(bedtools intersect -a "$tmp_merged_bed_file" -b "$tmp_merged_bed_file" -wo |\
   mlr --tsv --skip-comments --implicit-csv-header --headerless-csv-output filter '!($2 == $8 && $3 == $9 && $6 == $12)'\
   then count)
number_of_intersections=$((twice_number_of_intersections/2))
echo "\"${prefix}number_of_intersections_ignoring_strand_and_ignoring_self\": $number_of_intersections,"

# shellcheck disable=SC1010
closest_distance=\
$(bedtools closest -io -d -s -t first -a "$tmp_bed_file" -b "$tmp_bed_file" |\
    awk '{if($NF > 0){print $NF}}' |\
    mlr --itsv --ojson --implicit-csv-header rename 1,region_length\
      then stats1 -a mean,stddev,count,min,max,p10,p30,p50,p70,p90 -f region_length)

# shellcheck disable=SC1010
closest_distance_ignoring_strand=\
$(bedtools closest -io -d -t first -a "$tmp_bed_file" -b "$tmp_bed_file" |\
    awk '{if($NF > 0){print $NF}}' |\
    mlr --itsv --ojson --implicit-csv-header rename 1,region_length\
      then stats1 -a mean,stddev,count,min,max,p10,p30,p50,p70,p90 -f region_length)

echo "\"${prefix}closest_distance_statistics_per_region_same_strand_bp\": $closest_distance,"
echo "\"${prefix}closest_distance_statistics_per_region_ignoring_strand_bp\": $closest_distance_ignoring_strand,"
echo "\"comment\": \"NOTE: closest distance counts may be lower than overall region count \
because there may be no neighbours for some intervals\","
echo "\"${prefix}number_of_regions_merge_same_strand\": $(wc -l < "$tmp_merged_bed_file"),"
echo "\"${prefix}number_of_+_regions_merge_same_strand\": $(< "$tmp_merged_bed_file" awk 'BEGIN{a=0}$6 == "+"{a+=1}END{print a}'),"
echo "\"${prefix}number_of_-_regions_merge_same_strand\": $(< "$tmp_merged_bed_file" awk 'BEGIN{a=0}$6 == "-"{a+=1}END{print a}'),"

# if fasta file is available, count the number of bases in the bed file merged with same strand
if [ -f "$fasta_file" ]; then
  tmp_base_count_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  python count_bases_windows_given_fasta_and_bed.py -i "$tmp_merged_bed_file" -g "$fasta_file" |\
    awk -v OFS='\t' 'BEGIN{a=0;c=0;g=0;t=0}{a+=$((NF-4));c+=$((NF-3));g+=$((NF-2));t+=$((NF-1))}END{print a, c, g, t}' \
    > "$tmp_base_count_file"

  echo "\"${prefix}total_A_merge_same_strand\": $(< "$tmp_base_count_file" awk '{print $1}'),"
  echo "\"${prefix}total_C_merge_same_strand\": $(< "$tmp_base_count_file" awk '{print $2}'),"
  echo "\"${prefix}total_G_merge_same_strand\": $(< "$tmp_base_count_file" awk '{print $3}'),"
  echo "\"${prefix}total_T_merge_same_strand\": $(< "$tmp_base_count_file" awk '{print $4}'),"
  echo "\"${prefix}ratio_of_AT_to_total_merge_same_strand\": $(< "$tmp_base_count_file" awk '{print ($1 + $4)/($1 + $2 + $3 + $4)}'),"

  rm "$tmp_base_count_file"

fi

echo "\"${prefix}total_length_bp_merge_same_strand\": $(< "$tmp_merged_bed_file" awk 'BEGIN{a=0}{a+=$3-$2}END{print a}'),"
echo "\"${prefix}total_length_bp_merge_ignore_strand\": $(< "$tmp_merged_bed_file_ignore_strand" awk 'BEGIN{a=0}{a+=$3-$2}END{print a}'),"

# now, calculate some statistics about columns with numeric values
# shellcheck disable=SC2016
# shellcheck disable=SC1010
echo "\"column_value_statistics\":" \
  "$(grep -E -v '^#|^browser|^track' "$bed_file" |\
      mlr --itsv --ojson --implicit-csv-header put '$length=$3-$2' \
         then stats1 -a mean,median,stddev,min,max -f 5,7,8,9,10,length  | jq -r '.[]'),"

echo "\"comment_stats\": " \
  "\"In the calculation above, we've assumed columns 5,7,8,9,10 exist and are numeric. " \
  "If this assumption is not correct for any field, then ignore the corresponding output above\""

echo "}"

# remove temporary files
rm "$tmp_bed_file" "$tmp_merged_bed_file" "$tmp_merged_bed_file_ignore_strand"
if [ "$bed_file" == "$tmp_bed_file_stdin" ]; then
  rm "$tmp_bed_file_stdin"
fi
