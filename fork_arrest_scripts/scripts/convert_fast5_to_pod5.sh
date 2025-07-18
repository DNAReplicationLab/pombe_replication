#!/bin/bash
#SBATCH --mem=120G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J fast5toPod5
#SBATCH --mail-type=END,FAIL
#SBATCH --time=48:00:00

# goal
# ====
# convert fast5 files to pod5 files

# usage
# =====
# bash convert_fast5_to_pod5.sh <fast5_dir> <output_pod5>
# Can also use sbatch
# We will look recursively for fast5 files in the directory specified

set -e

# load pod5 package
source load_package.sh -pod5

# load input arguments
fast5_dir=${1:-}
output_pod5=${2:-}

print_usage() {
  >&2 echo "$1"
  >&2 echo "Usage: bash convert_fast5_to_pod5.sh <fast5_dir> <output_pod5>"
  >&2 echo "Can also use sbatch"
  >&2 echo "We will look recursively for fast5 files in the directory specified"
  exit 1
}

# check if the input directory exists
if [ ! -d "$fast5_dir" ]; then
  print_usage "fast5_dir does not exist"
fi

# check that the output file name is not empty
if [ -z "$output_pod5" ]; then
    print_usage "Error: output_pod5 is empty"
fi

# perform conversion
pod5 convert from_fast5 "$fast5_dir" -o "$output_pod5" --recursive --force-overwrite