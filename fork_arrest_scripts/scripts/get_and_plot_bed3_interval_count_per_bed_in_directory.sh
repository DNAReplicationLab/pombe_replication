#!/bin/bash

# goal
# -----
# Given a bed file (in BED3 format, call this A) and a directory of bed files, count and plot number of intervals in A
# that overlap with each bed file in the directory.

# usage and logic
#----------------
# bash get_and_plot_bed3_interval_count_per_bed_in_directory.sh <input_dir> <bed_file> <output_file> <json_options>
# input_dir - directory with bed files
# bed_file - bed file in BED3 format
# output_file - two output files are produced, a text file named output_file and a png file named output_file.png
# json_options - optional json file with options for plotting, can be omitted, see the appropriate section below for
#                details.

# output file description
# -----------------------
# The text output file is a tab separated file with two columns, bed_file and interval_count.
# bed_file - name of the bed file in the input directory
# interval_count - number of intervals in the input bed file that overlap with the bed file in the input directory
# The image output file is a scatter plot of interval counts per bed file.
# Some log messages are printed to the standard output.

# json file description
# ---------------------
# The json file looks like the following. One or both fields can be omitted.
# {
#   "bed_file_regex": "*.bed",
#   "bed_sort_file_name_numerical": false
# }
# bed_file_regex - regex pattern to match bed files in the input directory, defaults to *.bed means match all bed files
#                  if not provided.
# bed_sort_file_name_numerical - sort bed files in the input directory numerically, defaults to false if not provided.

# stop execution if any command fails
set -e

# load python
source load_package.sh -bedtools -R -miller -jq -python

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# print log message
echo "Receiving and checking command line arguments..."

# get command line arguments
input_dir=${1:-}
bed_file=${2:-}
output_file=${3:-}
json_options=${4:-}

# check if bed file exists
if [ ! -f "$bed_file" ]; then
  echo "Bed file does not exist: $bed_file"
  exit 1
fi

# check if bed file is in the correct format
if [ ! "$(< "$bed_file" python validate_bed_format.py --max-three-columns)" == "valid" ]; then
    >&2 echo "Error: input file bed_file is not in the correct bed format."
    exit 1;
fi

# check if json options file exists
is_json_options=false
if [ -f "$json_options" ]; then
  is_json_options=true
fi

# if available, then set some options
if [ "$is_json_options" == true ]; then
  bed_file_regex=$(< "$json_options" jq -r '.bed_file_regex // "*.bed"')
  bed_sort_file_name_numerical=$(< "$json_options" jq -r '.bed_sort_file_name_numerical // false')
else
  bed_file_regex="*.bed"
  bed_sort_file_name_numerical=false
fi

# check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: input directory does not exist"
    exit 1
fi

# make the parent directory of the output file if it does not exist
output_dir=$(dirname "$output_file")
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# begin processing
# ----------------
# print log message
echo "Processing..."

# prepare input bed file
tmp_bed_file=$(mktemp -p "${config[scratchDir]:-}"/tmp XXXXXX)
< "$bed_file" python print_valid_data_bed_lines.py |\
  awk -v OFS="\t" '{if($2!=$3){print $1, $2, $3}}' |\
    sort -k 1,1 -k2,2n |\
      bedtools merge -i stdin  > "$tmp_bed_file"

# get list of bed files in the input directory
if [ "$bed_sort_file_name_numerical" = true ]; then
  bed_file_list=$(find "$input_dir" -name "$bed_file_regex" | sort -V)
else
  bed_file_list=$(find "$input_dir" -name "$bed_file_regex" | sort)
fi

# print log message
echo "Begin counting intervals..."

# write header to output file
echo -e "bed_file\tinterval_count" > "$output_file"

# perform calculation
for current_bed_file in $bed_file_list; do

  if [ -f "$current_bed_file" ]; then
    # get the name of the bed file
    current_bed_file_name=$(basename "$current_bed_file")

    # count number of intervals in the bed file that overlap with the input bed file
    < "$current_bed_file" python print_valid_data_bed_lines.py |\
        awk -v OFS="\t" '{if($2!=$3){print $1, $2, $3}}' |\
          sort -k 1,1 -k2,2n |\
            bedtools merge -i stdin |\
              bedtools intersect -a "$tmp_bed_file" -b - -wa -u | wc -l |\
                awk -v OFS='\t' -v bed_file="$current_bed_file_name" '{print bed_file, $1}' >> "$output_file"

  else
    break;
  fi

done;

# check that there are at least 2 lines in the output file, or throw an error
if [ "$(wc -l < "$output_file")" -lt 2 ]; then
  >&2 echo "Error: output file has less than 2 lines."
  rm "$tmp_bed_file"
  rm "$output_file"
  exit 1;
fi

echo "Finished counting intervals."
echo "Plotting..."

# plot counts per bed file
# NOTE: we are repurposing the mean column for the count as the plotting script expects a column named mean
# shellcheck disable=SC1010
< "$output_file" mlr --tsv --skip-comments rename bed_file,x_label,interval_count,mean |\
  Rscript plotting_and_short_analyses/plot_scatter_categorical_x.R \
    "$output_file".png "Bed file name" "Count"

# remove temporary bed file
rm "$tmp_bed_file"

# print log message
echo "Finished plotting."