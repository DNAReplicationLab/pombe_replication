#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotRead
#SBATCH --mail-type=NONE
#SBATCH --mail-user=
#SBATCH --time 0:39:59
#SBATCH --constraint=""

# introduction
# --------

# program plots reads with raw data (required), windowed data (optional),
# model curve(s) (optional), and fork feature annotation (optional).

# NOTE: the script has several features, so please ensure you read the comments before executing the script.

# plotting data
# -------------
# a sample program execution line and its meaning follows

# bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir 300 0.5 1100 fasta 0.7_0.1_10000_300_10000_20000
# can also use sbatch instead of bash.

# the line above extracts data from the read with readid = readID from sample.mod.bam.
# Data runs on chrII from 1000 to 20000 and is windowed in 300 thymidines after thresholding
# i.e. calling T as BrdU if probability > 0.5. Windows are such that 1100 on the reference genome
# falls on a window boundary.
# The plot is sent to output_dir and has the name plot_readID.png.
# plot data is also sent to output_dir.
# fasta is the reference genome.
# Model curve is the best-fit sigmoid with a high, low of 0.7, 0.1, midpoint at 10000 and width of
# 300, that runs from 10000 to 20000.
# Multiple model curves can be plotted by simply appending
# more strings in the format "high_low_offset_width_modelStart_modelEnd" (see section below).

# variations on the command above to add or subtract plot features
# ----------------------------------------------------------------

# to get just the raw data, use
# bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir

# to get just raw and windowed data, use
# bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir 300 0.5 1100

# to plot paused forks, simply pass two "high_low_offset_width_modelStart_modelEnd" strings
# separated by a space i.e.
# bash plot_read.sh ... 0.7_0.1_15000_300_10000_20000 0.7_0.1_5000_300_20000_30000
# this will plot one curve running from x=10000 to x=20000, and another curve from x=20000 to x=30000.
# the discontinuity in the transition at x=20000 is a model-identified pause.
# NOTE: in the string above, '...' does not literally mean ... It means there are other parameters that are not shown.

# to plot fork features, pass a fork feature annotation file after the last parameter i.e.
# bash plot_read.sh ... fasta 0.7_0.1_15000_300_10000_20000 ... fork_feature_annotation_file
# bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir fork_feature_annotation_file
# and similar usages.
# Such an annotation is optional.
# fork_feature_annotation_file is a three-column, space-separated file without headers.
# the columns are start, end, label.
# start and end refer to positions on the reference genome.
# label can be 'origin', 'leftFork', 'rightFork', 'termination' or 'pause'.
# start must always be less than or equal to end i.e. don't switch coordinates for left/right fork.

# if you want to plot multiple forks or multiple paused forks, read this section carefully.
# to show multiple forks, we need to introduce breaks in the model curve to signify that we
# are showing separate forks. In R, that is achieved using NA values supplied to geom_path.
# so, we are going to use a special format of parameter string to tell our program to introduce
# breaks in the model curve. The format is "12000_NA", which means insert an NA at x = 12000.
# sample usage:
# bash plot_read.sh ... 0.7_0.1_10000_300_10000_20000 25000_NA 0.7_0.1_35000_300_30000_40000
# NOTE: in the string above, '...' does not literally mean ... It means there are other parameters that are not shown.
# in the usage above, a model curve is plotted from 10000 to 20000, and a model curve from 30000 to 40000,
# and these curves are spatially separated i.e. there's no black line running from 20000 to 30000 because
# an NA has been inserted at 25000 (any number from 20000 to 30000 will do, no need to use the midpoint).

# stop execution if any command fails
set -e

# set output directory, making it if it doesn't exist
mkdir -p "$6";
op_dir=$(cd "$6"; pwd)

# load R, python
pwd=$(pwd)
config_dir=..
cd "$config_dir"
source load_package.sh -R -python

# initialize annotation file name with a random string
annotation_file=$(openssl rand -hex 16)

# write function to check if parameter strings are in the valid format
is_valid_parameter_string() {
    local input_string="$1"
    IFS='_' read -ra entries <<< "$input_string"
    local is_valid=1

    for entry in "${entries[@]}"; do
        if [[ ! "$entry" =~ ^-?[0-9]+(\.[0-9]+)?$ ]] && [[ "$entry" != "NA" ]]; then
            is_valid=0
            break
        fi
    done

    if [[ $is_valid -eq 1 ]]; then
        echo "valid"
    else
        echo "invalid"
    fi
}

# set filenames
data_file="$op_dir"/plot_data_"$2"
plot_file="$op_dir"/plot_"$2".png
mod_bam_file="$op_dir"/plot_"$2".bam

# gather data
{

  # get information about the read
  bash get_information_from_read_id.sh -o "$mod_bam_file" "$1" "$2";

	# get raw data
	python get_raw_data_from_modBAM.py "$mod_bam_file" "$2" "$3" "$4" "$5" |\
	  awk '{print $1 " " $2 " " $2+1 " " $3 " rawDetect"}';

  # window it
  if [[ $# -gt 6 ]] && [[ $7 =~ ^[0-9]+$ ]]
  then
    python get_raw_data_from_modBAM.py "$mod_bam_file" "$2" "$3" "$4" "$5" |\
      sed '1i\detectIndex\tposOnRef\tval' |\
      python get_mean_brdU_window.py --window "$7" --thres "$8" --forceWinBoundaryAtPos "$9" |\
      awk '{print $1 " " $3 " " $4 " " $2 " winDetect"}';
	fi

	# add model curve(s)
	for i in "${@:11}"
	do

	  if [[ $(is_valid_parameter_string "$i") == "invalid" ]]; then
      continue
    fi

	  check_if_NA=$(echo "$i" | cut -d "_" -f 2)

	  if [ "$check_if_NA" == "NA" ]
	  then
	    position_of_NA=$(echo "$i" | cut -d "_" -f 1)
	    echo "$2 $position_of_NA $(echo "$position_of_NA" + 1 | bc) NA model";
	  else
      start_fork=$(echo "$i" | cut -d "_" -f 5)
      end_fork=$(echo "$i" | cut -d "_" -f 6)
      param_str=$(echo "$i" | cut -d "_" -f 1,2,3,4)
      python get_model_values_at_raw_data_coords.py sigmoidal "$3" "$start_fork" "$end_fork" "$param_str" "${10}" |\
              awk '{print $1 " " $1+1 " " $2 " model"}' |\
              sed "s/^/$2 /";
    fi
  done

} > "$data_file"

cd "$pwd"
# get the last argument on the command line
# treat that as the annotation file, if such a file exists and annotation file is not already set
for last in "$@"; do true; done

if [ -f "$last" ] && [ ! -f "$annotation_file" ]
then
	annotation_file="$last";
fi

# plot data
if [ -e "$annotation_file" ]
then
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$data_file" "$plot_file" "$annotation_file"
else
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$data_file" "$plot_file"
fi