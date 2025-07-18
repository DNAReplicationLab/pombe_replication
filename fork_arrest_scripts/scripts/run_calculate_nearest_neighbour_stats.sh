#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J runThPauseDistanceStatsAndPlot
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# This script calculates the nearest neighbour distance between pauses and plots the distribution of these distances,
# along with the theoretical distribution of distances between pauses.

# usage
#------
# bash run_calculate_nearest_neighbour_stats.sh $pause_file $sensitivity_file $fasta_fai $op_dir [$max_distance]
# NOTE: [] denotes an optional argument. To provide it, remove the square brackets and type the value.
# $pause_file: The pause file in our usual format
# $sensitivity_file: The sensitivity file in our usual format
# $fasta_fai: The fasta index file
# $op_dir: The output directory
# $max_distance: The maximum distance to consider between pauses. Default is 10000.

# outputs
# -------
# A table with the nearest neighbour distances between pauses and the theoretical distribution of
# distances between pauses and the corresponding plot is saved in the output directory.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -bedtools -miller -R

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 something"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "The suffix _optional means that the parameter is optional."
    exit 1;
fi

# assign arguments to variables
pause_file=$1
sensitivity_file=$2
fasta_fai=$3
op_dir=$4
max_distance=${5:-10000}
bin_size=100

# check that the first three arguments are files and exist
for file in "$pause_file" "$sensitivity_file" "$fasta_fai"; do
  if [ ! -f "$file" ]; then
    >&2 echo "ERROR: $file is not a file or does not exist."
    exit 1
  fi
done

# check that the input pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# check that max distance is a multiple of bin size
if [ $((max_distance % bin_size)) -ne 0 ]; then
  >&2 echo "ERROR: max_distance must be a multiple of $bin_size, but $max_distance is not a multiple of $bin_size."
  exit 1;
fi
# check that max distance is greater than bin size
# shellcheck disable=SC2086
if [ $max_distance -le $bin_size ]; then
  >&2 echo "ERROR: max_distance must be greater than $bin_size, but $max_distance is not greater than $bin_size."
  exit 1;
fi

# make output directory if it does not exist
mkdir -p "$op_dir"
op_dir=$(realpath "$op_dir")

# calculate the theoretical distribution of distances between pauses
bash run_calculate_theoretical_pause_neighbour_distance_from_sensitivity.sh "$sensitivity_file" \
  "$max_distance" "$fasta_fai" "$bin_size" > "$tmpDir"/theoretical_pause_neighbour_distance.txt

# convert pause file to bed format
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
  grep -v '^#' |  sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,$4}' > "$tmpDir"/pause.bed

# ensure that the lines in the file above are unique or raise an error
if [ "$(sort -u "$tmpDir"/pause.bed | wc -l)" -ne "$(wc -l < "$tmpDir"/pause.bed)" ]; then
  >&2 echo "ERROR: pause file contains duplicate entries."
  exit 1
fi

# use bedtools closest to calculate the distance between each pause and its nearest neighbour
bedtools closest -N -a "$tmpDir"/pause.bed -b "$tmpDir"/pause.bed -d -t first |\
 awk -F"\t" '$NF >= 0 {print $NF}' |\
  mlr --tsv --implicit-csv-header histogram -f 1 --lo 0 --hi "$max_distance" \
    --nbins $((max_distance/bin_size)) > "$tmpDir"/closest_pause.bed

# join the two files
# shellcheck disable=SC1010
# shellcheck disable=SC2016
mlr --itsv --ocsv --skip-comments put '$bin_lo=fmtnum($bin_lower_limit,"%d")' \
  then join -f "$tmpDir"/closest_pause.bed -j bin_lo "$tmpDir"/theoretical_pause_neighbour_distance.txt |\
 insert_calling_script_header "$@" > "$op_dir"/combined_dist_table.csv

# plot the result
cd plotting_and_short_analyses
Rscript plot_nearest_neighbour_histogram.R "$op_dir"/combined_dist_table.csv "$op_dir"/combined_dist_table.png

# remove temporary directory
rm -r "$tmpDir"