#!/bin/bash

# goal
# -----
# Subset four bam files to retain only the reads that are present in the read_ids file

# usage
#------
# bash subset_four_bam_files.sh $bam $mod_bam $mod_bam_left $mod_bam_right $read_ids $op_dir
# NOTE: if the bam files are large, this process can take a long time and/or be memory intensive.
#       We leave it to the user to set the appropriate resources for the job using slurm or otherwise.
# bam: path to the bam file, set to /dev/null if not available/unused
# mod_bam: path to a mod bam file, set to /dev/null if not available/unused
# mod_bam_left: path to a mod bam left file, set to /dev/null if not available/unused
# mod_bam_right: path to a mod bam right file, set to /dev/null if not available/unused
# read_ids: path to the read ids file, a plain text file with one read id per line and no header or comments
# op_dir: path to the output directory, will be created if it does not exist

# outputs
# -------
# A json object with the paths to the subset bam files is written to stdout.
# Four files (and associated .bai indices) are created in the output directory called
# bam_subset.bam, mod_bam_subset.bam, mod_bam_left_subset.bam, mod_bam_right_subset.bam
# If any file is set to /dev/null in the input, it is not created, and the output json contains the path /dev/null
# If any such named files already exist, they are overwritten.

# stop execution if any command fails
set -e

# load packages
source load_package.sh -samtools

# assign arguments to variables
bam_file="$1"
mod_bam="$2"
mod_bam_left="$3"
mod_bam_right="$4"
read_id_processed="$5"
op_dir="$6"

# check that the correct number of arguments were provided
if [ "$#" -lt 6 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash subset_four_bam_files.sh bam mod_bam mod_bam_left mod_bam_right read_ids op_dir"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# check that the read_ids file exists
if [ ! -f "$read_id_processed" ]; then
    >&2 echo "ERROR: read_ids file not found"
    exit 1;
fi

# make the output directory if it does not exist
mkdir -p "$op_dir"

# Function to create subsets of bam files
create_bam_subset() {
  local bam_file="$1"
  local subset_prefix="$2"
  local bam_subset="/dev/null"

  if [[ -f "$bam_file" ]]; then
    bam_subset="${op_dir}/${subset_prefix}.bam"
    samtools view -b -N "$read_id_processed" "$bam_file" > "${bam_subset}.tmp"
    samtools sort -o "$bam_subset" "${bam_subset}.tmp"
    samtools index "$bam_subset"
    rm "${bam_subset}.tmp"
  fi

  echo "$bam_subset"
}

# * Make subsets of bam files wherever possible
# ** Create temporary files to capture the output of background tasks
tmp_file_1=$(mktemp -p "$op_dir")
tmp_file_2=$(mktemp -p "$op_dir")
tmp_file_3=$(mktemp -p "$op_dir")
tmp_file_4=$(mktemp -p "$op_dir")

# ** Make subsets of bam files wherever possible and capture their output in temporary files
(create_bam_subset "$bam_file" "bam_subset" > "$tmp_file_1") &
(create_bam_subset "$mod_bam" "mod_bam_subset" > "$tmp_file_2") &
(create_bam_subset "$mod_bam_left" "mod_bam_left_subset" > "$tmp_file_3") &
(create_bam_subset "$mod_bam_right" "mod_bam_right_subset" > "$tmp_file_4") &

# ** Wait for all background tasks to finish
wait

# Output a json object with the paths to the subset bam files
echo "{"
echo "\"bam_subset\": \"$(cat "$tmp_file_1")\","
echo "\"mod_bam_subset\": \"$(cat "$tmp_file_2")\","
echo "\"mod_bam_left_subset\": \"$(cat "$tmp_file_3")\","
echo "\"mod_bam_right_subset\": \"$(cat "$tmp_file_4")\""
echo "}"

# clean up temporary files
rm "$tmp_file_1" "$tmp_file_2" "$tmp_file_3" "$tmp_file_4"