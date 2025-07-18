#!/bin/bash

# goal
# -----
# we add some columns to the pauseFile from rDNA analysis so that it works with our plotting routines

# usage
#------
# bash rDNA_prepare_pauseFile_for_plotting.sh input_file output_file
# names are self-explanatory.
# need not use sbatch as this is a very quick script

# outputs
# -------
# output sent to output_file with some columns added from input_file

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -python

# load configuration
source config.sh

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
input_file=${1:-}
output_file=${2:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
  echo >&2 "ERROR: incorrect number of arguments provided"
  echo >&2 "usage: bash run_rDNA_prepare_pauseFile_for_plotting.sh input_file output_file"
  echo >&2 "For more details on what these parameters mean, see the comments in the script."
  exit 1
fi

# check that input file exists
if [ ! -f "$input_file" ]; then
  echo >&2 "ERROR: input file does not exist"
  echo >&2 "input_file: $input_file"
  exit 1
fi

# check that it is a valid pause file
if [ ! "$(< "$input_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# shellcheck disable=SC2016
< "$input_file" mlr --itsv --otsv --skip-comments \
  put '$start = $boundary_left_of_step;
       $end = $boundary_right_of_step;
       $paramStrLeft = $dens_at_left_of_step . "_" . $dens_at_left_of_step . "_" . "NA_NA";
       $paramStrRight = $dens_at_right_of_step . "_" . $dens_at_right_of_step . "_" . "NA_NA";' |\
  insert_calling_script_header "$@" > "$output_file"
