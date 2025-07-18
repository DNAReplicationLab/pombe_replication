#!/bin/bash
#SBATCH --mem=20G
#SBATCH -p ei-medium
#SBATCH -J calcAutoCorr
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# goal
# ----
# Calculate autocorrelation per read of modBAM reads, then isolate reads with correlation lengths > window size,
# then plot some of these reads grouping by correlation length.
# Our overall goal is to see if correlations can be used to isolate reads with analogue density gradients in them.
# (Our experimental protocols produce many kinds of reads and only those with gradients are of interest.)

# usage
# -----
# sbatch run_calculate_modBAM_autocorrelation.sh $mod_bam $output_dir
# where:
# $mod_bam = path to modBAM file
# $output_dir = path to output directory, preferably an empty directory, will be created if it doesn't exist
# can use bash if a cluster is not available or if this is a short job

# output
# ------
# Many plots and text files are sent to the output directory. These include:
# - files showing autocorrelation function vs lag window size per chunk per read\
# - files showing autocorrelation length per chunk per read
# - files showing max autocorrelation length per read
# - folders with plots associated with reads grouped by max correlation length

# fail on error
set -e

# load configuration
source config.sh

# set main and plotting directories
main_dir=$(pwd)
plotting_scripts_dir="$main_dir"/plotting_and_short_analyses

# load python
source load_package.sh -python -samtools -miller

# load input arguments
mod_bam=${1:-/dev/null}
output_dir=$2

# check that at least 2 arguments were provided
if [ $# -lt 2 ]; then
    echo "Usage: sbatch $0 mod_bam output_dir"
    echo "For more information, see the header of this script"
    exit 1
fi

# check that mod_bam exists
if [ ! -f "$mod_bam" ] || [ ! -f "$mod_bam".bai ]; then
    echo "mod_bam and/or index not found: $mod_bam"
    exit 1
fi

# make output directory if necessary
output_dir=$(realpath "$output_dir")
mkdir -p "$output_dir"

# calculate autocorrelation
echo "Calculating autocorrelation"
samtools view "$mod_bam" |\
  python calculate_modBAM_autocorrelation.py --chunk 8000 --chunkOverlap 4000 \
    --textFilePrefix "$output_dir"/measure > "$output_dir"/measure_autocorr.txt

# histogram of max autocorrelation lengths
echo "Calculating histogram of max autocorrelation lengths"
mlr --tsv histogram -f max_auto_correlation_length --lo 0 --hi 10 --nbins 10 "$output_dir"/measure_summary_read.txt \
  > "$output_dir"/max_auto_correlation_length_histogram.tsv

# make plots for reads with max autocorrelation lengths of 1 to 7
echo "Making plots for reads grouping by max autocorrelation lengths"
for i in $(seq 1 7); do
  # make a directory for this max autocorrelation length
  plot_dir="$output_dir"/plot_max_auto_correlation_length_"$i"
  mkdir -p "$plot_dir"

  # filter reads with max autocorrelation length of i and select up to 10 reads at random
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  mlr --tsv --skip-comments --headerless-csv-output filter '$max_auto_correlation_length == '"$i" \
    then cut -f detectIndex "$output_dir"/measure_summary_read.txt | shuf | head -n 10 |\
    awk -F"_" '{print $1}' > "$plot_dir"/read_ids.txt

  # now plot these reads
  cd "$plotting_scripts_dir"
  blank=/dev/null
  sbatch plot_n_reads.sh "$plot_dir"/read_ids.txt "$mod_bam" "$blank" "$blank" "$blank" "$blank" \
     "$blank" "$plot_dir"
  cd "$main_dir"

done