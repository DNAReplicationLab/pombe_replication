#!/bin/bash

# goal
# -----
# Given a directory of bed files, extract a numeric column in each bed file and report its mean and sd per bed file.

# usage
#------
# bash get_and_plot_mean_bed_value_per_bed_in_directory.sh <input_dir> <column_number> \
#   <output_dir> <json_options>
# can use sbatch instead of bash if script takes a while to run. if so, then set the sbatch options suitably
# NOTE: \ at the end of the line is for line continuation. It can be removed if the command is written in one line.
# NOTE: [] means optional.
# <input_dir> - directory with bed files
# <column_number> - column number to calculate mean and sd. 1-based. default is 5 (5th column i.e. score)
# <output_dir> - directory to save output
# <json_options> - (default unused) optional json file with options for plotting. Format is explained below.

# format of json_options
# ----------------------
# the file should be a json file with the following format (fields can be in any order):
# {
#    "bed_file_regex": "bed_LS_*.bed",
#    "bed_sort_file_name_numerical": true,
#    "output_file_prefix": "blah_blah",
#    "y_label_plot": "blah_blah"
# }
# <bed_file_regex> - regex to match bed files in <input_dir>, default is *.bed
# <bed_sort_file_name_numerical> - if true, sort bed files in <input_dir> numerically by file name, otherwise sort
#                                  alphabetically. default is false.
# <output_file_prefix> - prefix for output files. default is "mean_signal_per_bed"
# <y_label_plot> - y-axis label for plot. default is "Median replication time (min)" (due to historical reasons)

# outputs
# -------
# Plot(s), and a table with mean, sd, and count of signal per bed file are sent to <output_dir>

# stop execution if any command fails
set -e

# load packages
source load_package.sh -jq -R -miller -python

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# get command line arguments
input_dir=${1:-}
column_number=${2:-5}
output_dir=${3:-}
json_options=${4:-}

# check if number of arguments is correct
if [ "$#" -lt 3 ]; then
    >&2 echo "Error: incorrect number of arguments"
    >&2 echo "usage: bash get_and_plot_mean_bed_value_per_bed_in_directory.sh <input_dir> <column_number> <output_dir> <json_options>"
    >&2 echo "<input_dir> - directory with bed files"
    >&2 echo "<column_number> - column number to calculate mean and sd. 1-based. default is 5 (5th column i.e. score)"
    >&2 echo "<output_dir> - directory to save output"
    >&2 echo "<json_options> - optional json file with options for plotting, can be omitted. See script for details."
    exit 1
fi

# check if input directory exists
if [ ! -d "$input_dir" ]; then
    >&2 echo "Error: input directory does not exist"
    exit 1
fi

# check if output directory exists, if not create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    output_dir=$(realpath "$output_dir")
fi

# check that column number is a number
if ! [[ "$column_number" =~ ^[0-9]+$ ]]; then
    >&2 echo "Error: column number is not a number"
    exit 1
fi

# check that column number is larger than 3
if [ "$column_number" -le 3 ]; then
    >&2 echo "Error: column number must be larger than 3"
    exit 1
fi

# check if json options file exists
is_json_options=false
if [ -f "$json_options" ]; then
  is_json_options=true
fi

# if available, then set some options
if [ "$is_json_options" == true ]; then

  bed_file_regex=$(< "$json_options" jq -r '.bed_file_regex')
  bed_sort_file_name_numerical=$(< "$json_options" jq -r '.bed_sort_file_name_numerical')
  output_file_prefix=$(< "$json_options" jq -r '.output_file_prefix')
  y_label_plot=$(< "$json_options" jq -r '.y_label_plot')

  if [ "$bed_file_regex" == null ]; then
    bed_file_regex="*.bed"
  fi

  if [ "$bed_sort_file_name_numerical" == null ]; then
    bed_sort_file_name_numerical=false
  fi

  if [ "$output_file_prefix" == null ]; then
    output_file_prefix="mean_signal_per_bed"
  fi

  if [ "$y_label_plot" == null ]; then
    y_label_plot="Median replication time (min)"
  fi

else
  bed_file_regex="*.bed"
  bed_sort_file_name_numerical=false
  output_file_prefix="mean_signal_per_bed"
  y_label_plot="Median replication time (min)"
fi

# loop over bed files
header_line=1
{

  echo "# bed_file_regex: $bed_file_regex"
  echo "# bed_sort_file_name_numerical: $bed_sort_file_name_numerical"
  echo "# output_file_prefix: $output_file_prefix"

  if [ "$bed_sort_file_name_numerical" = true ]; then
      bed_file_list=$(find "$input_dir" -name "$bed_file_regex" | sort -V)
    else
      bed_file_list=$(find "$input_dir" -name "$bed_file_regex" | sort)
  fi

  for bed_file in $bed_file_list; do

    if [ -f "$bed_file" ]; then

      # perform calculations
      tmp_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
      # shellcheck disable=SC1010
      < "$bed_file" python print_valid_data_bed_lines.py | awk '$2 < $3' |\
        mlr --itsv --ojson --implicit-csv-header rename "$column_number",signal \
          then stats1 -a mean,stddev,count -f signal | grep -v -E "^\[|^\]" > "$tmp_file"

      mean=$(< "$tmp_file" jq -r '.signal_mean')
      sd=$(< "$tmp_file" jq -r '.signal_stddev')
      count=$(< "$tmp_file" jq -r '.signal_count')

      # if mean, sd or count are non-numeric, then set them to NA
      if [ -z "$mean" ] || [[ ! "$mean" =~ ^[0-9.e+-]+$ ]]; then
        mean="NA"
      fi
      if [ -z "$sd" ] || [[ ! "$sd" =~ ^[0-9.e+-]+$ ]]; then
        sd="NA"
      fi
      if [ -z "$count" ] || [[ ! "$count" =~ ^[0-9]+$ ]]; then
        count="NA"
      fi

      if [ "$header_line" -eq 1 ]; then
        echo -e "x_label\tmean\tsd\tcount"
        header_line=0
      fi

      bed_file_name=$(basename "$bed_file")

      echo -e "$bed_file_name\t$mean\t$sd\t$count"
      rm "$tmp_file";

    else
      break;
    fi
  done;

} > "$output_dir"/"$output_file_prefix".txt

# plot the data
< "$output_dir"/"$output_file_prefix".txt Rscript plotting_and_short_analyses/plot_error_bar_categorical_x.R \
  "$output_dir"/"$output_file_prefix".png "Bed file name" "$y_label_plot"

# plot the standard error of the mean (SEM)
# NOTE: we are repurposing the sd column for the SEM as the plotting script expects a column named sd
# shellcheck disable=SC1010,SC2016
< "$output_dir"/"$output_file_prefix".txt mlr --tsv --skip-comments \
  put -e 'if($count > 1){$y = $sd/sqrt($count)}else{$y = "NA"}' \
  then cut -f x_label,mean,y then rename y,sd |\
  Rscript plotting_and_short_analyses/plot_error_bar_categorical_x.R \
    "$output_dir"/"$output_file_prefix".sem.png "Bed file name" "$y_label_plot"

# plot counts per bed file
# NOTE: we are repurposing the mean column for the count as the plotting script expects a column named mean
# shellcheck disable=SC1010
< "$output_dir"/"$output_file_prefix".txt mlr --tsv --skip-comments cut -f x_label,count then rename count,mean |\
  Rscript plotting_and_short_analyses/plot_scatter_categorical_x.R \
    "$output_dir"/"$output_file_prefix".count.png "Bed file name" "Count"

# plot sums per bed file
# NOTE: we are repurposing the mean column for the sum as the plotting script expects a column named mean
# shellcheck disable=SC1010
# shellcheck disable=SC2016
< "$output_dir"/"$output_file_prefix".txt mlr --tsv --skip-comments cut -f x_label,count,mean \
  then rename mean,mn then put -e 'if($mn == "NA" || $count == "NA"){$sum = "NA"} else {$sum = $count * $mn}'\
  then rename sum,mean |\
    Rscript plotting_and_short_analyses/plot_scatter_categorical_x.R \
      "$output_dir"/"$output_file_prefix".sum.png "Bed file name" "Sum"