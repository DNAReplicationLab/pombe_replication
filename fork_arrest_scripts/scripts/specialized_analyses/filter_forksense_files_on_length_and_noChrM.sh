#!/bin/bash

# goal
# ----
# filter out forksense features that are too short or lie on chrM

# usage
# -----
# bash <program.sh> $forksense_dir $min_fork_length $min_align_length $output_dir $chrM_name
# don't bother using sbatch, just run it directly as the script is very fast
# $forksense_dir is the directory containing the forksense files with standard names and in standard format
# $min_fork_length (optional) is the minimum fork length in kb above which forks are retained, default is 10
# $min_align_length (optional) is the minimum alignment length in kb above which forks on that alignment are retained,
#     default 30
# $output_dir (optional) is the directory to which the filtered forksense files are written.
#     default is $forksense_dir
# $chrM_name (optional) is the name of the mitochondrial chromosome, default is chrM

# script fails if any command fails
set -e

# set input parameters
forksense_dir=$1
min_fork_length=${2:-10}
min_align_length=${3:-30}
output_dir=${4:-$forksense_dir}
chrM_name=${5:-chrM}

# check at least one argument was provided
if [ $# -lt 1 ]
then
    >&2 echo "Usage: bash <program.sh> <forksense_dir> <min_fork_length> <min_align_length> <output_dir> <chrM_name>"
    >&2 echo "min_fork_length (in kb) and min_align_length (in kb) are optional and default to 10 and 30 respectively."
    >&2 echo "output_dir is optional and defaults to forksense_dir"
    >&2 echo "chrM_name is optional and defaults to chrM"
    exit 1;
fi

# check left and right fork files exist
left_fork_file="$forksense_dir"/leftForks_DNAscent_forkSense.bed
right_fork_file="$forksense_dir"/rightForks_DNAscent_forkSense.bed
if [[ ! -f "$left_fork_file" ]] || [[ ! -f "$right_fork_file" ]]; then
  >&2 echo "left fork file or right fork file do not exist!"
  exit 1;
fi

# check these files have eight columns
if [[ $(awk '{print NF}' "$left_fork_file" | sort -nu) -ne 8 ]] ||\
   [[ $(awk '{print NF}' "$right_fork_file" | sort -nu) -ne 8 ]]; then
  >&2 echo "left fork file or right fork file do not have eight columns!"
  exit 1;
fi

# capitalize the first letter of chrM_name and store it in a variable
chrM_name_capitalized=$(echo "$chrM_name" | awk '{print toupper(substr($0,1,1)) substr($0,2)}')

# set left and right fork files which contain a subset of the forks
mkdir -p "$output_dir"
left_fork_subset_file="$output_dir"/\
leftForks_DNAscent_forkSense.forkLen.${min_fork_length}kb.alignLen.${min_align_length}kb.no"$chrM_name_capitalized".bed
right_fork_subset_file="$output_dir"/\
rightForks_DNAscent_forkSense.forkLen.${min_fork_length}kb.alignLen.${min_align_length}kb.no"$chrM_name_capitalized".bed

# multiply min_fork_length and min_align_length by 1000 to get them in bp
min_fork_length_bp=$((min_fork_length*1000))
min_align_length_bp=$((min_align_length*1000))

# perform filtering on left and right forks
awk -v min_fork_length="$min_fork_length_bp" -v min_align_length="$min_align_length_bp" \
    'BEGIN{OFS="\t"}{if ($3 - $2 >= min_fork_length && $7 - $6 >= min_align_length && $1 != "'"$chrM_name"'") print $0}' \
    "$left_fork_file" > "$left_fork_subset_file"

awk -v min_fork_length="$min_fork_length_bp" -v min_align_length="$min_align_length_bp" \
    'BEGIN{OFS="\t"}{if ($3 - $2 >= min_fork_length && $7 - $6 >= min_align_length && $1 != "'"$chrM_name"'") print $0}' \
    "$right_fork_file" > "$right_fork_subset_file"

# if origin file exists, filter it too
origin_file="$forksense_dir"/origins_DNAscent_forkSense.bed
if ! [[ -f "$origin_file" ]]; then
  exit 0;
fi

# load the package miller and configuration variables
cd ..
source load_package.sh -miller > /dev/null
source config.sh

# set temporary directory and make it
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# join the left and right forks to the origins
cd plotting_and_short_analyses;
bash join_dnascent_origins_and_forks.sh "$origin_file" "$left_fork_subset_file" "$right_fork_subset_file" "$tmpDir" \
  "$min_fork_length" "$min_align_length" > /dev/null

# process the joined file so that we can output the origins in the format used by forksense
origin_subset_file="$output_dir"/\
origins_DNAscent_forkSense.forkLen.${min_fork_length}kb.alignLen.${min_align_length}kb.no"$chrM_name_capitalized".bed

# shellcheck disable=SC2016
# shellcheck disable=SC1010
< "$tmpDir"/combine_origin_leftForks_rightForks.tsv > "$origin_subset_file" \
  mlr --itsv --ocsv --ofs ' ' --skip-comments --headerless-csv-output \
  put 'b = splitax($detectIndex,"_");
       $contig = b[2];
       $st = b[3];
       $en = b[4];
       $strand = b[5];'\
  then cut -o -f contig,OL2,OL3,read_id,contig,st,en,strand\
  then filter '$contig != "'"$chrM_name"'"'

# remove the temporary file
rm "$tmpDir"/combine_origin_leftForks_rightForks.tsv;

# output paths to the filtered files in json format
echo "{"
echo "  \"left_fork_subset_file\": \"$left_fork_subset_file\","
echo "  \"right_fork_subset_file\": \"$right_fork_subset_file\","
echo "  \"origin_subset_file\": \"$origin_subset_file\""
echo "}"