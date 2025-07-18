#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J runIntExactBed
#SBATCH --mail-type=END,FAIL
#SBATCH --time=2:59:59
#SBATCH --constraint=""

# goal
# -----
# intersect two bed files in an exact manner
# for more details, refer to the documentation of the python script

# usage
#------
# sbatch run_intersect_exact_bed_files.sh <bed1> <bed2> <output_file>
# can use bash instead of sbatch but the script might take a few minutes to run

# inputs
# ------
# bed1: first bed file
# bed2: second bed file

# outputs
# -------
# output_file: output bed file

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -python

# load configuration
source config.sh

# ensure three arguments are provided
if [ "$#" -ne 3 ]; then
    >&2 echo "Illegal number of parameters"
    >&2 echo "Usage: sbatch run_intersect_exact_bed_files.sh <bed1> <bed2> <output_file>"
    >&2 echo "bed1: first bed file"
    >&2 echo "bed2: second bed file"
    >&2 echo "output_file: output bed file"
    exit 1;
fi

# check that the first and second bed files exist
if [ ! -f "$1" ]; then
    >&2 echo "Input file $1 does not exist"
    exit 1;
fi

if [ ! -f "$2" ]; then
    >&2 echo "Input file $2 does not exist"
    exit 1;
fi

{
  # print some info about the run
  echo "# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>";

  # print the command used to run the script
  echo "# command: $0 $*";

  python intersect_exact_bed_files.py "$1" "$2";

} > "$3"
