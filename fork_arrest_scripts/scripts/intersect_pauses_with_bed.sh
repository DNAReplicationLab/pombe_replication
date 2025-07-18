#!/bin/bash

# Goal
# ----
# Given a set of genomic intervals and pause-calling results, group calls based on intersection with the intervals.

# Details
# -------
# Basically, the three categories of pause-calling results are:
# - pauseSite lies in region
# - pauseSite not in region but fork lies in region
# - no pause detected but fork lies in region
# We want to produce these three files and their union in a fourth file.
# All pause calls are associated with a fork and it is assumed that a fork can have zero or one pause call.

# Usage
# -----
# Example usages:
# bash intersect_pauses_with_bed.sh -s $region_file $pause_file $all_forks $no_pause $elsewhere_pause $region_pause
# bash intersect_pauses_with_bed.sh -S $region_file $pause_file $all_forks $no_pause $elsewhere_pause $region_pause
# bash intersect_pauses_with_bed.sh -f $region_file $pause_file $all_forks $no_pause $elsewhere_pause $region_pause
# bash intersect_pauses_with_bed.sh -fS $region_file $pause_file $all_forks $no_pause $elsewhere_pause $region_pause
#   region_file: bed file with regions of interest.
#   pause_file: pause-calling results, tab-separated with column names.
#               Must have columns pauseSite pauseDuration keep* detectIndex
#                 - pauseSite is a location on the reference genome in bp. (contig can be inferred from detectIndex)
#                 - pauseDuration is the duration of the pause in bp.
#                 - keep* are columns with boolean values indicating whether the pause
#                   passed that particular criterion, set to 'True' or 'False'
#                 - detectIndex is a string in the format
#                   readid_contig_start_end_orientation_direction_startFork_endFork.
#                   It means a fork with direction L or R was detected at position startFork to endFork on the
#                   alignment given by readid_contig_start_end_orientation.
#                   orientation can be fwd/rev.
#                   start < end and startFork < endFork irrespective of orientation or direction.
#   all_forks: output file, see the section Details.
#   no_pause: output file, see the section Details.
#   elsewhere_pause: output file, see the section Details.
#   region_pause: output file, see the section Details.
#   -s option: if used, intersect only if same strand.
#   -S option: if used, intersect only if opposite strand.
#   -f option: if used, use fork direction L/R to set strand orientation -/+. Use with -s or -S, not by itself.

# Output
# ------
# See the section Details.

# fail upon error
set -e

# Initialize variables for optional inputs
small_s=false
capital_s=false
small_f=false

# Parse optional arguments
while getopts "sSf" opt; do
  case $opt in
    s)
      small_s=true
      ;;
    S)
      capital_s=true
      ;;
    f)
      small_f=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Check if both -s and -S are set
if [ "$small_s" == "true" ] && [ "$capital_s" == "true" ]; then
  >&2 echo "Error: Both -s and -S options are set. Please use only one of them."
  exit 1
fi

# Check that if -f is set, then either -s or -S should be set as well
if [ "$small_f" == "true" ] && [ "$small_s" == "false" ] && [ "$capital_s" == "false" ]; then
  >&2 echo "Error: -f option is set but neither -s nor -S is set. Please use -f with either -s or -S."
  exit 1
fi

# Shift arguments to exclude parsed options
shift $((OPTIND-1))

# need 6 arguments
if [ $# -lt 6 ]; then
  >&2 echo "Error: incorrect number of arguments"
  >&2 echo "Usage: bash $0 [-s] [-S] [-f] \$region_file \$pause_file \$all_forks \$no_pause \$elsewhere_pause \$region_pause"
  >&2 echo "  [-s] [-S] [-f] means use nothing or a combination of -s, -S, -f , -fS, -f -S, etc."
  >&2 echo "  -s option: if used, intersect only if same strand."
  >&2 echo "  -S option: if used, intersect only if opposite strand."
  >&2 echo "  -f option: if used, use fork direction L/R as the strand orientation -/+."
  >&2 echo "  NOTE: to use strand information, the bed file must have a 6th column with strand information."
  >&2 echo "  NOTE: do not use -f by itself. Use with -s or -S."
  >&2 echo "  NOTE: do not use both -s and -S."
  >&2 echo "  region_file: bed file with regions of interest."
  >&2 echo "  pause_file: pause-calling results, tab-separated with column names."
  >&2 echo "  all_forks, no_pause, elsewhere_pause, region_pause: output files."
  exit 1;
fi

region_file=$1
pause_file=$2
op_all_pause_file=$3 # here, "all" means pauses whether we have confidence in them or not
op_no_pause_file=$4 # forks with no pauses in them
op_elsewhere_pause_file=$5 # forks that pass through region but pause elsewhere
op_region_pause_file=$6 # forks that pass through region and pause in it

# check that input files exist
if [ ! -f "$region_file" ] || [ ! -f "$pause_file" ]; then
  >&2 echo "Error: one or both input files do not exist."
  exit 1;
fi

# load user information
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# load bedtools, python, and miller
source load_package.sh -bedtools -python -miller

# load git repo labels
source load_git_repo_labels.sh

# depending on strand specific usages, check that the bed file has the right number of columns
if [ "$small_s" == "true" ] || [ "$capital_s" == "true" ] || [ "$small_f" == "true" ]; then
  if [ ! "$(< "$region_file" python validate_bed_format.py --no-dot-strand --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed file must have at least six valid columns with sixth column = +/- when using -s/-S/-f."
    exit 1;
  fi
else
  if [ ! "$(< "$region_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed file must have at least three valid columns."
    exit 1;
  fi
fi

# check that the pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# function to extract header and column names
function extract_header_and_column_names() {
  local file="$1"
  head -200 "$file" | grep '^#' | cat
  # CAUTION: the above will fail if there are more than 200 lines of comments, but what is the likelihood of that?
  grep -m 1 -v '^#' "$file"
}

# set flags for strand-specificity
strand_specific_option=""
if [ "$small_s" == "true" ]; then
  strand_specific_option="-s"
elif [ "$capital_s" == "true" ]; then
  strand_specific_option="-S"
fi

# set flags for fork direction
fork_direction_option=""
if [ "$small_f" == "true" ]; then
  fork_direction_option="--LRtoPlusMinus"
fi

# create a file with all the pause information associated with forks that overlap with the regions of interest
# ============================================================================================================

# get fork coordinates in bed format
tmp_file_forks=$(mktemp -p "${config[scratchDir]:-}"/tmp)
< "$pause_file" python convert_pause_file_to_bed.py --outputForks $fork_direction_option > "$tmp_file_forks"

# intersect tmp_file_forks with region_file
{
  extract_header_and_column_names "$pause_file";

  bedtools intersect -u -wa $strand_specific_option -a "$tmp_file_forks" -b "$region_file" | cut -f7-

} > "$op_all_pause_file"

# create a file with all information about forks with no pauses that overlap with the regions of interest
# =======================================================================================================

{
  extract_header_and_column_names "$op_all_pause_file"

  # get those forks that have no pause
  < "$op_all_pause_file" python convert_pause_file_to_bed.py |\
       grep -v '^#' | awk '{if($5 == 0){print $0}}' | cut -f7-

} > "$op_no_pause_file"

# create a file with all the pause information associated with pauses that lie in the region of interest
# =======================================================================================================

# get pause coordinates in bed format, keeping only those with confidence 1000
tmp_file_pauses=$(mktemp -p "${config[scratchDir]:-}"/tmp)
< "$op_all_pause_file" python convert_pause_file_to_bed.py $fork_direction_option  |\
     awk '{if($5 == 1000){print $0}}' > "$tmp_file_pauses"

# intersect tmp_file_pauses with region_file
{
  extract_header_and_column_names "$op_all_pause_file"

  # finally, perform the intersect
  bedtools intersect -u -wa $strand_specific_option -a "$tmp_file_pauses" -b "$region_file" | cut -f7-

} > "$op_region_pause_file"

# create a file with all the pause information associated with pauses that lie outside the region of interest
# ===========================================================================================================

# invert the intersection between tmp_file_pauses and region_file
{
  extract_header_and_column_names "$op_all_pause_file"

  # finally, perform the invert-intersect
  bedtools intersect -v -wa $strand_specific_option -a "$tmp_file_pauses" -b "$region_file" | cut -f7-

} > "$op_elsewhere_pause_file"

# delete temporary files
rm "$tmp_file_forks" "$tmp_file_pauses"