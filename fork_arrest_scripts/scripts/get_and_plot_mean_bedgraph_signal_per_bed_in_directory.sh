#!/bin/bash

# goal
# -----
# Given a bedgraph with a signal and directory of bed files, calculate and plot the mean and sd of signal per bed

# usage
#------
# bash get_and_plot_mean_bedgraph_signal_per_bed_in_directory.sh <input_dir> <bedgraph_file> \
#   <output_dir> <json_options> <intersect_options>
# can use sbatch instead of bash if script takes a while to run. if so, then set the sbatch options suitably
# NOTE: \ at the end of the line is for line continuation. It can be removed if the command is written in one line.
# NOTE: [] means optional.
# <input_dir> - directory with bed files
# <bedgraph_file> - bedgraph file with signal in the 4th column. This is a space-separated file with no column names,
#                   comments starting with #, and columns: contig, start, end, signal.
# <output_dir> - directory to save output
# <json_options> - (default unused) optional json file with options for plotting. Format is explained below.
# <intersect_options> - (default is "..i") string with options for calculate_bed_overlap_bedgraph_signal_stats.sh.
#                       string must have three characters.
#                       the first can be p or m or . (for plus or minus or no strand for the bedgraph),
#                       the second can be s or S or . (for same or opposite or no strand for intersections with bed),
#                       and the third can be i or e (for intensive or extensive overlap).
#                       See calculate_bed_overlap_bedgraph_signal_stats.sh for details.
#                       e.g.: pSi, mse, psi, mSi, ..e etc.

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
source load_package.sh -jq -R -miller

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# get command line arguments
input_dir=${1:-}
bedgraph_file=${2:-}
output_dir=${3:-}
json_options=${4:-}
intersect_options=${5:-..i}

# check if number of arguments is correct
if [ "$#" -lt 3 ]; then
    echo "Error: incorrect number of arguments"
    echo "usage: bash get_mean_bedgraph_signal_per_bed_in_directory.sh <input_dir> <bedgraph_file> <output_dir> <json_options> <intersect_options>"
    echo "<input_dir> - directory with bed files"
    echo "<bedgraph_file> - bedgraph file with signal in the 4th column."
    echo "<output_dir> - directory to save output"
    echo "<json_options> - optional json file with options for plotting, can be omitted. See script for details."
    echo "<intersect_options> - optional options for strand and bedgraph value, can be omitted. See header of this script for details."
    exit 1
fi

# check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: input directory does not exist"
    exit 1
fi

# check if bedgraph file exists
if [ ! -f "$bedgraph_file" ]; then
    echo "Error: bedgraph file does not exist"
    exit 1
fi

# check if output directory exists, if not create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    output_dir=$(realpath "$output_dir")
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

# parse intersect_options. Remove any periods, and add a hyphen before every character and store in a new array
intersect_options=${intersect_options//./}
intersect_options=${intersect_options//e/} # we are removing extensive overlap, as it is the default option
                                           # for calculate_bed_overlap_bedgraph_signal_stats.sh
intersect_options_array=()
for (( i=0; i<${#intersect_options}; i++ )); do
  intersect_options_array+=("-${intersect_options:$i:1}")
done

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
      bash calculate_bed_overlap_bedgraph_signal_stats.sh "${intersect_options_array[@]}" "$bed_file"\
        "$bedgraph_file" > "$tmp_file"
      mean=$(< "$tmp_file" jq -r '.overlap_mean')
      sd=$(< "$tmp_file" jq -r '.overlap_stddev')
      count=$(< "$tmp_file" jq -r '.overlap_count')

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