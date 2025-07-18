#!/usr/bin/env bash

# goal
# =====
# Given a bed file, and a bedgraph with a signal value, calculate stats of signal for all overlaps with bed file.
# By default, we assume the signal is extensive e.g. if a bed file entry overlaps 50% of the length of a bedgraph entry,
# we take 50% of the signal for that instance of overlap.
# The user can specify that the signal is intensive e.g. if a bed file entry overlaps 50% of the length of a bedgraph
# entry, we take 100% of the signal for that instance of overlap.
# An example of an extensive quantity is number of features like pauses in an interval.
# So, if we consider half the interval size for example, we should expect to pick up half the number of pauses.
# An example of an intensive quantity is median replication time of an interval.
# So, if we consider half the interval size for example, we should expect to pick up the same median replication time,
# and a calculation that halves the median replication time is incorrect.
# Extensive and intensive are terms from physics: https://en.wikipedia.org/wiki/Intensive_and_extensive_properties

# usage
# =====
# bash calculate_bed_overlap_bedgraph_signal_stats.sh [-s|-S] [-p|-m] [-i] <bed file> <bedgraph file>
# NOTE: <> indicates required arguments which must be provided in the order shown.
#       Please fill in the path to the files you want to use in place of the placeholders (e.g. <bed file>).
# NOTE: [] indicates optional arguments.
#       -s: if used, intersect only if same strand.
#       -S: if used, intersect only if opposite strand.
#       -p: if used, interpret the bedgraph as having data corresponding to the plus strand.
#       -m: if used, interpret the bedgraph as having data corresponding to the minus strand.
#       -i: if used, interpret the bedgraph as having intensive data. By default, interpret data as extensive.
#       To do a stranded intersection, use -s or -S and -p or -m.
#       Otherwise, we perform the intersect ignoring the strand.
# NOTE: to use strand information, the bed file must have a 6th column with strand information.
#       If strand options are set and the bed file does not have strand information, the script will exit with an error.
# bed file: bed file of interest
# bedgraph file: bedgraph file of interest, columns are contig, start, end, signal without column names,
#                and space-separated.

# some examples
# =============
# bash calculate_bed_overlap_bedgraph_signal_stats.sh /path/to/bed/file /path/to/bedgraph/file
#   - ignore strand information and get stats of signal for intersections
# bash calculate_bed_overlap_bedgraph_signal_stats.sh -s -p /path/to/bed/file /path/to/bedgraph/file
#   - only intersect if same strand and interpret bedgraph as having data corresponding to plus strand
# bash calculate_bed_overlap_bedgraph_signal_stats.sh -S -m /path/to/bed/file /path/to/bedgraph/file
#   - only intersect if opposite strand and interpret bedgraph as having data corresponding to minus strand

# output
# ======
# A json object is output to stdout and contains many fields that report statistics.
# {
#    "overlap_sum": 123.45,
#    "overlap_mean": 12.34,
#    "overlap_blah": 11,
#    "overlap_blah2": 22,
# }
# and so on.
# NOTE: blah and blah2 are just placeholders for other fields that are output.
# NOTE: and so on means that there are many other fields that are output.

# fail on error
set -Eeuo pipefail

# process optional arguments
small_s=false
capital_s=false
plus=false
minus=false
intensive=false

# keep track of number of flags set out of -s, -S, -p, -m
flag_count=0

# process flags
while getopts "sSpmi" opt; do
    case $opt in
        s)
            small_s=true
            flag_count=$((flag_count+1))
            ;;
        S)
            capital_s=true
            flag_count=$((flag_count+1))
            ;;
        p)
            plus=true
            flag_count=$((flag_count+1))
            ;;
        m)
            minus=true
            flag_count=$((flag_count+1))
            ;;
        i)
            intensive=true
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

# if both -p and -m are set, exit
if [ "$plus" == "true" ] && [ "$minus" == "true" ]; then
    >&2 echo "ERROR: both -p and -m options are set. Please use only one of them."
    exit 1
fi

# ensure that if -s/-S is set, then -p/-m is also set and vice versa
if [ "$flag_count" -eq 1 ] || [ "$flag_count" -eq 3 ] || [ "$flag_count" -eq 4 ]; then
    >&2 echo "ERROR: Set one of -s/-S and one of -p/-m or none of these flags at all."
    exit 1;
fi

# shift arguments to exclude parsed options
shift $((OPTIND-1))

# check that the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash calculate_bed_overlap_bedgraph_signal_stats.sh [-s|-S] [-p|-m] [-i] <bed file> <bedgraph file>"
    >&2 echo "NOTE: <> indicates required arguments which must be provided in the order shown."
    >&2 echo "NOTE: [] indicates optional arguments."
    >&2 echo "To do a stranded intersect use -s or -S to specify same or opposite strand, "
    >&2 echo "    and -p or -m to specify if the bedgraph has data corresponding to the plus or minus strand."
    >&2 echo "Otherwise, we perform the intersect ignoring the strand."
    >&2 echo "NOTE: to use strand information, the bed file must have a 6th column with strand information."
    >&2 echo "-i is an optional flag to interpret the bedgraph as having intensive data. By default, interpret data as extensive."
    >&2 echo "Extensive means a quantity that scales with interval size e.g. number of pauses."
    >&2 echo "Intensive means a quantity that does not scale with interval size e.g. median replication time."
    >&2 echo "For more explanations, refer to the script header."

    exit 1
fi

# get input arguments
bed_file=$1
bedgraph_file=$2

# check that these files exist
if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file does not exist"
    exit 1
fi

if [ ! -f "$bedgraph_file" ]; then
    >&2 echo "ERROR: bedgraph file does not exist"
    exit 1
fi

# load packages
source load_package.sh -bedtools -python -miller

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# depending on strand specific usages, check that the bed file has the right columns
if [ "$small_s" == "true" ] || [ "$capital_s" == "true" ]; then

  if [ ! "$(< "$bed_file" python validate_bed_format.py --no-dot-strand --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least six valid columns with sixth column = +/- when using -s/-S."
    exit 1;
  fi

else

  if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed files must have at least three valid columns."
    exit 1;
  fi

fi

# convert bedgraph to bed file depending on strand specific usages
tmp_bed_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)

if [ "$plus" == "true" ]; then
  strand_bedgraph="+"
elif [ "$minus" == "true" ]; then
  strand_bedgraph="-"
else
  strand_bedgraph="."
fi

grep -E -v "^#|^browser|^track" "$bedgraph_file" |\
  awk -v strand="$strand_bedgraph" 'BEGIN{OFS="\t"}{print $1, $2, $3, "blank", $4, strand}' > "$tmp_bed_file"

# check that the tmp bed file is of the right format
if [ ! "$(< "$tmp_bed_file" python validate_bed_format.py --exactly-six-columns --allow-float-score --never-zero-base)" == "valid" ]; then
  >&2 echo "Error: bedgraph is not of the correct format."
  exit 1;
fi

# perform calculation
tmp_calc_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
{
  if [ "$small_s" == "true" ] && [ "$flag_count" -eq 2 ]; then
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -s -i stdin -c 6 -o distinct | awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1, $2, $3, "blank", 1000, $4}}' |\
        bedtools intersect -s -a stdin -b "$tmp_bed_file" -wo
  elif [ "$capital_s" == "true" ] && [ "$flag_count" -eq 2 ]; then
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -s -i stdin -c 6 -o distinct | awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1, $2, $3, "blank", 1000, $4}}' |\
        bedtools intersect -S -a stdin -b "$tmp_bed_file" -wo
  else
    sort -k 1,1 -k2,2n "$bed_file" |\
      bedtools merge -i stdin | awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1, $2, $3, "blank", 1000, "."}}' |\
        bedtools intersect -a stdin -b "$tmp_bed_file" -wo
  fi
} > "$tmp_calc_file"

# if there are no lines in the calculation file, then print NA for all stats, except count which is 0
if [ ! -s "$tmp_calc_file" ]; then
  echo '{'
  for key in "mean" "stddev" "sum" "min" "max" "p10" "p30" "p50" "p70" "p90"; do
    echo "  \"overlap_$key\": \"NA\","
  done
  echo '  "overlap_count": 0 '
  echo '}'
  exit 0
fi

# shellcheck disable=SC1010
grep -v '^#' "$tmp_calc_file" |\
  awk -v intensive="$intensive" -v OFS="\t" '{if (intensive == "true") {print $1, $2, $3, $6, $13, $11 * $13} else {print $1, $2, $3, $6, 0, $11 * $13/($9 - $8)}}' |\
  bedtools groupby -g 1,2,3,4 -c 5,6 -o sum,sum  |\
  awk -v intensive="$intensive" '{if (intensive == "true") {print $6/$5} else {print $6}}' |\
  mlr --itsv --ojson --implicit-csv-header rename 1,overlap\
      then stats1 -a mean,stddev,sum,count,min,max,p10,p30,p50,p70,p90 -f overlap |\
  grep -v -E '^\[|^\]'

# remove temporary files
rm "$tmp_bed_file"
rm "$tmp_calc_file"