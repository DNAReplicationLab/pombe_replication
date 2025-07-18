#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotReadFromReadData
#SBATCH --mail-type=NONE
#SBATCH --mail-user=
#SBATCH --time 0:19:59
#SBATCH --constraint=""

# description
# -----------
# plot data from plot_data_* files in a directory
# usage: bash <script_name.sh> $directory
# can use sbatch in place of bash
# directory: path to a directory that contains a bunch of plot data files with names in the format plot_data_readid
#            where read id is the read id of a molecule from the sequencing experiment. The file contains tab-separated
#            data. For details on its format, refer to the comments in the plotting script
#            plot_one_read_w_win_or_model_if_needed.R.

# why is this script needed?
# --------------------------
# This script is used to re-run plot making, so why are we doing this?
# The answer is for publications, posters etc. one might want to change the plot options.
# So, we just dump the plot_data files in a folder, adjust options in the R script,
# and re-run this script to make the plots.
# We don't want to make publication-ready scripts in our pipeline as that would involve a lot of
# customisation-related parameters as input to the pipeline, which would make the pipeline harder to use.
# The workflow is: the pipeline quickly tells you if the analysis is worth pursuing, and then you can
# use this script to make the plots with the options you want.

# stop execution if any command fails
set -e

# load R, python
pwd=$(pwd)
config_dir=..
cd "$config_dir"
source load_package.sh -R
cd "$pwd"

# loop through all files in the directory and make plots suitably
op_dir=$1

for file in "$op_dir"/*
do
  # get basename of the file
  filename=$(basename "$file")

  # check if file starts with "plot_data_"
  if [[ $filename == plot_data_* ]]; then
    # cut "plot_data_" from the start of the filename
    suffix="${filename#plot_data_}"
    # plot data
    if [ -f "$op_dir"/features_"$suffix" ]; then
      Rscript ./plot_one_read_w_win_or_model_if_needed.R "$op_dir"/plot_data_"$suffix" "$op_dir"/plot_"$suffix".png \
        "$op_dir"/features_"$suffix";
    else
      Rscript ./plot_one_read_w_win_or_model_if_needed.R "$op_dir"/plot_data_"$suffix" "$op_dir"/plot_"$suffix".png;
    fi
  fi
done