#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 2
#SBATCH -p ei-medium
#SBATCH -J getMinMaxAnalogueSlidingWindow
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Given mod bam file and list of forks in the detectIndex format, print a min and max analogue density measured using
# sliding windows per fork

# usage
#------
# sbatch get_min_max_analogue_sliding_window.sh <detectIndex_file> <modbam> <n>
# detectIndex file should be tab-separated, can have comments starting with '#', must contain column names
#     as the first non-comment line, must contain a column called detectIndex in the format
#     readid_contig_start_end_strand_forkDirection_startFork_endFork.
#       - strand is either + or -
#       - forkDirection is either L or R
#       - start <= end irrespective of strand,
#       - startFork <= endFork irrespective of forkDirection,
#       - contig can contain underscores,
#       - startFork to endFork must lie within start to end.
# modbam file should be a modified bam file
# n (optional) is the number of thymidines per sliding window, defaults to 300

# outputs
# -------
# sent to standard output
# if you use sbatch -o something.txt ... (where ... means the rest of the command),
# then the output will be in something.txt

# stop execution if any command fails
set -e

# load packages
source load_package.sh -python -miller

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
  echo >&2 "ERROR: incorrect number of arguments provided"
  echo >&2 "usage: sbatch get_min_max_analogue_sliding_window.sh <detectIndex_file> <modbam> <n>"
  echo >&2 "detectIndex file should be tab-separated and have a detectIndex column. For further details on file format, see the script."
  echo >&2 "modbam file should be a modified bam file"
  echo >&2 "n (optional) is the number of thymidines per sliding window, defaults to 300"
  exit 1
fi

# check that the input files exist
if [ ! -f "$1" ]; then
  echo >&2 "ERROR: $1 does not exist"
  exit 1
fi

if [ ! -f "$2" ]; then
  echo >&2 "ERROR: $2 does not exist"
  exit 1
fi

# check that n is a positive integer
n=${3:-300}
if ! [[ "$n" =~ ^[0-9]+$ ]]; then
  echo >&2 "ERROR: n must be an integer >= 0"
  exit 1
fi

suffix="$n"T_sliding_window_along_fork

# make a temporary file
tmp_file=$(mktemp)

# split the detectIndex file into 50 files using tmp_file as a prefix
# shellcheck disable=SC1010
mlr --tsv --skip-comments rename detectIndex,alt_read_id then cat "$1" |\
  mlr --tsv split -m 3 --prefix "$tmp_file"

# do processing per split file
for splitFile in "$tmp_file"_*
do
  # shellcheck disable=SC1010
  < "$splitFile" \
    python get_raw_data_from_modBAM.py --piped-regions --fork-info-in-alt-read-id "$2" |\
      sed '1idetectIndex\tposOnRef\tval' |\
        python get_mean_brdU_window.py --thres 0.5 --sliding --window "$n" |\
          sed '1idetectIndex mean_val start end' |\
            mlr --icsv --ifs ' ' --otsv --skip-comments stats1 -a min,max -f mean_val -g detectIndex then\
              rename mean_val_min,min_"$suffix",mean_val_max,max_"$suffix" > "$splitFile"_processed &
done

# wait for calculations to finish
wait;

# concatenate results and print to standard output
mlr --tsv --skip-comments cat "$tmp_file"_*_processed

# remove temporary files
rm "$tmp_file"*