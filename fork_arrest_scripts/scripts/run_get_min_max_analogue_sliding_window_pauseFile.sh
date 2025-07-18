#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 20
#SBATCH -p ei-medium
#SBATCH -J getMinMaxAnalogueSlidingWinFromPauseFile
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Given a pause_file and a modbam file, add two columns of min and max mean analogue density using sliding windows

# usage
#------
# sbatch run_get_min_max_analogue_sliding_window_pauseFile.sh pause_file modbam_file
# - pause_file is a file containing pauses, it must pass checks in validate_pause_format.py, so refer to that file
#   for more details on the format
# - modbam_file is a bam file containing modification information.

# outputs
# -------
# - sent to stdout.
# - if you want it sent to a file, use sbatch -o something.txt ... where something.txt is the name of the output file
#     and ... is the rest of the sbatch command

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -miller

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


# check that the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch run_get_min_max_analogue_sliding_window_pauseFile.sh pause_file modbam_file"
    >&2 echo "pause_file is a file containing pauses, it must pass checks in validate_pause_format.py, so refer to that file for format information."
    >&2 echo "modbam_file is a bam file containing modification information."
    exit 1;
fi

# assign arguments to variables
pause_file="$1"
modbam_file="$2"

# check if these files exist
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause_file does not exist"
    exit 1;
fi

if [ ! -f "$modbam_file" ]; then
    >&2 echo "ERROR: modbam_file does not exist"
    exit 1;
fi

# check that the pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# create two temporary file
tmp_file="$(mktemp)"
tmp_file_original="$(mktemp)"

# get the two new columns
bash get_min_max_analogue_sliding_window.sh "$pause_file" "$modbam_file" 300 | grep -v "^#" > "$tmp_file"

# make a pause file without comments
grep -v "^#" "$pause_file" > "$tmp_file_original"

# join the pause file with the temp file
mlr --itsv --otsv join -j detectIndex -f "$tmp_file" "$tmp_file_original" | insert_calling_script_header "$pause_file" "$modbam_file"

# remove the temporary files
rm "$tmp_file" "$tmp_file_original"