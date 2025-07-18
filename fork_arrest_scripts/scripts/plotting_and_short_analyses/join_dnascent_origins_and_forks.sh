#!/bin/bash

# goal
# -----
# DNAscent forkSense outputs four files: left fork, right fork, origins, terminations
# Each origin was called because from the site emerged a left and a right fork in a read.
# Now, we want to associate each origin call with the left and right fork that allowed such a call.

# usage
# -----
# bash join_dnascent_origins_and_forks.sh <origins> <left_forks> <right_forks> <output_dir>
#     <min_fork_length> <min_align_length>
# no need to use sbatch; the script is fast enough by itself.
# <origins> is the output of DNAscent forkSense, containing the origins
# <left_forks> is the output of DNAscent forkSense, containing the left forks
# <right_forks> is the output of DNAscent forkSense, containing the right forks
# <output_dir> is the directory where the output will be saved
# <min_fork_length> (optional) is the minimum fork length in kb above which forks are retained, default 10
# <min_align_length> (optional) is the minimum alignment length in kb above which forks on that alignment are retained,
#     default 30

# * the output file is tab separated with headers. it has the columns:
#   detectIndex, detectIndex_L, detectIndex_R, read_id, OL2, OL3, L2, L3, OR2, OR3, R2, R3
# * the interpretation of the columns is that a left fork b/w L2, L3 and a right fork b/w R2, R3 emerged from the origin
#   given by OL2, OL3. (OR2 and OR3 should be equal to OL2 and OL3 respectively).
# * L2 < L3, R2 < R3 i.e. the coordinates are sorted. In other words, the forks run from L3 -> L2 and R2 -> R3.
# * a sample entry is shown below, with each element displayed here on a separate line to help readability
#   although in the actual file, they'll all be on one line.
# * this is not real data, just an example
# 25e41ad9-89de-4698-a6be-eb12dd49a03d_chrII_10_55598_fwd
# 25e41ad9-89de-4698-a6be-eb12dd49a03d_chrII_10_55598_fwd_L_237_18500
# 25e41ad9-89de-4698-a6be-eb12dd49a03d_chrII_10_55598_fwd_R_19039_38672
# 25e41ad9-89de-4698-a6be-eb12dd49a03d
# 18500
# 19039
# 237
# 18500
# 18500
# 19039
# 19039
# 38672

# fail if any command fails
set -e

# set configuration directory
config_dir=..

# load packages, configuration variables, git variables
pwd=$(pwd)
cd "$config_dir" || exit;
source load_package.sh -miller
source config.sh;
source load_git_repo_labels.sh;
cd "$pwd" || exit;

# if the config variable doesn't exist, then exit
if [ -z "${config[scratchDir]}" ] || [ -z "${config[name]}" ]; then
    >&2 echo "Missing configuration value(s)";
    exit 1;
fi

# if the git details don't exist, then exit
if [ -z "${COMMITSTR}" ]; then
  >&2 echo "Missing git information";
  exit 1;
fi

# check that there are the correct number of arguments
if [ "$#" -lt 4 ]; then
  >&2 echo "Usage: bash join_dnascent_origins_and_forks.sh <origins> <left_forks> <right_forks> <output_dir> <min_fork_length> <min_align_length>"
  >&2 echo "  <origins> is the output of DNAscent forkSense, containing the origins"
  >&2 echo "  <left_forks> is the output of DNAscent forkSense, containing the left forks"
  >&2 echo "  <right_forks> is the output of DNAscent forkSense, containing the right forks"
  >&2 echo "  <output_dir> is the directory where the output will be saved"
  >&2 echo "  <min_fork_length> (optional) is the minimum fork length in kb above which forks are retained, default 10"
  >&2 echo "  <min_align_length> (optional) is the minimum alignment length in kb above which forks on that alignment are retained, default 30"
  exit 1;
fi

# receive arguments
origins=$1
left_forks=$2
right_forks=$3
op_dir=$4
min_fork_length=${5:-10}
min_align_length=${6:-30}

# check that the input files exist
if [ ! -f "$origins" ]; then
  >&2 echo "Origins file $origins does not exist"
  exit 1;
fi

if [ ! -f "$left_forks" ]; then
  >&2 echo "Left forks file $left_forks does not exist"
  exit 1;
fi

if [ ! -f "$right_forks" ]; then
  >&2 echo "Right forks file $right_forks does not exist"
  exit 1;
fi

# check that min_fork_length and min_align_length are integers greater than or equal to zero
if ! [[ "$min_fork_length" =~ ^[0-9]+$ ]] || ! [[ "$min_align_length" =~ ^[0-9]+$ ]]; then
  >&2 echo "min_fork_length and min_align_length must be integers greater than or equal to zero"
  exit 1;
fi

# multiply fork and alignment lengths by 1000 to convert to bp
min_fork_length_bp=$((min_fork_length * 1000))
min_align_length_bp=$((min_align_length * 1000))

# set temporary directory and files
# make temporary directory if need be
tmp_dir="$op_dir"/temp"$(openssl rand -hex 4)"
mkdir -p "$tmp_dir"

temp1="$tmp_dir"/temp1
temp2="$tmp_dir"/temp2
op_file="$op_dir"/combine_origin_leftForks_rightForks.tsv

# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --icsv --ifs ' ' --otsv --implicit-csv-header join -j 1,4,5,6,7,8 --lp OL --rp L -f "$origins" \
    then filter 'b = ($OL2 + $OL3)/2 - $L3; b <= 1000 && b >= 0' "$left_forks" > "$temp1"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --icsv --ifs ' ' --otsv --implicit-csv-header join -j 1,4,5,6,7,8 --lp OR --rp R -f "$origins" \
    then filter 'b = $R2 - ($OR2 + $OR3)/2; b <= 1000 && b >= 0' "$right_forks" > "$temp2"

{
  echo "# generated from ${COMMITSTR} at $(date) by ${config[name]}"
  echo "# input file: $origins"
  echo "# input file: $left_forks"
  echo "# input file: $right_forks"
  echo "# min_fork_length: $min_fork_length_bp"
  echo "# min_align_length: $min_align_length_bp"

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  mlr --itsv --otsv join -j 1,4,5,6,7,8 -f "$temp1" \
    then filter 'b = $R2 - $L3; b <= 2000 && b >= 0 &&
      $OL2 == $OR2 && $OL3 == $OR3 &&
      $7-$6 >= '$min_align_length_bp' &&
      $L3-$L2 >= '$min_fork_length_bp' &&
      $R3-$R2 >= '$min_fork_length_bp \
    then put '$detectIndex = $4. "_" . $5 . "_" . $6 . "_" . $7 . "_" . $8' \
    then put '$detectIndex_L = $detectIndex . "_L_" . $L2 . "_" . $L3' \
    then put '$detectIndex_R = $detectIndex . "_R_" . $R2 . "_" . $R3' \
    then cut -o -f detectIndex,detectIndex_L,detectIndex_R,4,OL2,OL3,L2,L3,OR2,OR3,R2,R3 \
    then rename 4,read_id \
    "$temp2";

}> "$op_file"

# remove temporary files
rm "$temp1" "$temp2";
rm -r "$tmp_dir";

# shellcheck disable=SC1010
echo "Histogram of number of origins per alignment (equivalent to per read id if one alignment per read id)"
# shellcheck disable=SC1010
mlr --itsv --opprint --skip-comments count -g detectIndex\
  then sort -n count then stats1 -a count -f detectIndex -g count "$op_file"