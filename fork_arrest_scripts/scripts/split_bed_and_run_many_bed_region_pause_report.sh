#!/bin/bash

# goal
# -----
# Given a bed file, split into many pieces by a column value and along each bed entry, then get pause report for
# each piece

# usage and inputs
#-----------------
# bash split_bed_and_run_many_bed_region_pause_report.sh input.json
# - no need of sbatch as we only perform simple steps before launching many slurm jobs here.
# - input.json is a json file with fields required by the two main scripts called here:
#   split_bed_file_by_column_value_and_along_length.sh and run_many_bed_region_pause_report.sh,
#   and an optional log_file field and a random_string_job_name field.
#   Please see the comments in those scripts for more details on what fields are required.
#   Example: let's say the first script requires the fields "a" and "b" and the second script requires the fields
#   "c" and "d". Then the input.json file should look like this (with fields in any order):
#   {
#     "a": "value1",
#     "b": "value2",
#     "c": "value3",
#     "d": "value4",
#     "log_file": "log.txt",
#     "random_string_job_name": "job_name"
#   }
#   log_file: (default /dev/null which means no log file) is the path to the log file to use.
#   random_string_job_name: (default "" i.e. unused) is a random string to use for the pause-report job name

# outputs
# -------
# This script just calls other scripts, so please look at those scripts for what their outputs are.
# Some messages are written to a log file and to the standard output i.e. on screen or to an SBATCH output file
# depending on how the job is run.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -jq -python

# load configuration
source config.sh

# print node information
bash print_node_information.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 input.json"
    >&2 echo "For more details on what parameters are needed, please see the comments in the script"
    exit 1;
fi

# load input parameters
input_json=$1

# check that the input file exists
if [ ! -f "$input_json" ]; then
    >&2 echo "ERROR: input file $input_json does not exist"
    exit 1;
fi

# make two temporary files
temp_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
temp_file_input="$temp_file".input.json

# process the input json file to expand any variables
< "$input_json" python expand_variables_in_json.py > "$temp_file_input"

# set the log file
log_file=$(jq -r '.log_file // "/dev/null"' "$temp_file_input")

# split the bed file
{
  echo "Splitting bed file according to the input json file"
  cat "$temp_file_input"
} | tee -a "$log_file"
bash split_bed_file_by_column_value_and_along_length.sh "$temp_file_input" > "$temp_file"

# combine input json with the temp file
jq --argfile temp_file "$temp_file" '."regions_of_interest" = $temp_file' "$temp_file_input" > "$temp_file".json

# call run_many_bed_region_pause_report.sh on each split bed file
{
  echo "Running many bed region pause report according to the input json file"
  cat "$temp_file".json
  bash run_many_bed_region_pause_report.sh "$temp_file".json
} | tee -a "$log_file"

# remove the temp file
rm "$temp_file" "$temp_file".json "$temp_file_input"