#!/bin/bash

# goal
# -----
# bring together statistics from several pause reports that are related e.g.: co-directional, head-on pauses for tRNA,
#    genes with different transcription rates, biological replicates etc.

# usage
#------
# bash collate_pause_report_statistics.sh input.json
# cat something.json | bash collate_pause_report_statistics.sh stdin
# cat something.json | bash collate_pause_report_statistics.sh -
# input.json: a json file containing a list of directories within which lie one or several pause report files

# inputs
# ------
# input.json must have the following format:
# {
#     "directories": ["/path/containing/pauses/report",
#                     "/another/path/containing/pauses/report",
#                     ...],
# }
# - ... means that you can add as many directories as you want
# - You can use any format for the input json file, as long as it is valid json.
# - The directories field contains a list of directories within which lie one or several pause report files,
#   which have filenames in the format "*_pauses_report.json" where * means a string of any length >= 0.

# outputs
# -------
# A many-column tab-separated table in plain text to standard output with column names such as "value", "desc"
# No comments are output

# stop execution if any command fails
set -e

# load packages
source load_package.sh -jq -miller

# load configuration variables
source config.sh

# set temporary directory and make it
tmp_dir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmp_dir"

# assign arguments to variables
input_json=${1:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash collate_pause_report_statistics.sh input.json"
    >&2 echo "usage: cat something.json | bash collate_pause_report_statistics.sh stdin"
    >&2 echo "usage: cat something.json | bash collate_pause_report_statistics.sh -"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# if input_json is stdin or -, then read from stdin, store in a temporary file, and set input_json to the temporary file
if [ "$input_json" == "stdin" ] || [ "$input_json" == "-" ]; then
  tmp_json=$(mktemp -p "$tmp_dir" tmp_json_XXXXXX.json)
  cat > "$tmp_json"
  input_json="$tmp_json"
fi

# check that the input json file exists
if [ ! -f "$input_json" ]; then
    >&2 echo "ERROR: input json file does not exist"
    exit 1;
fi

# extract the directories from the input json file, and append "*_pauses_report.json" to each directory
files=$(jq -r '.directories | .[]' "$input_json" | xargs -I {} echo {}/*_pauses_report.json)

# loop over the files, and extract the value and desc fields from each file
# then filter out any rows where value is null
# shellcheck disable=SC2016
{
  for file in $files; do
    cat "$file"
  done
} | jq '.[] | {value: .value, desc: .desc, dataset: .dataset,
               feature: .feature, division: .division,
               analysis_label: .analysis_label,
               relative_direction: .relative_direction}' | mlr --ijson --otsv filter '$value != null'

# clean up
rm -rf "$tmp_dir"