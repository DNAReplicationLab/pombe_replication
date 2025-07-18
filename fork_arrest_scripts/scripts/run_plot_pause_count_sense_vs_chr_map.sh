#!/bin/bash
#SBATCH --mem-per-cpu=100G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotPauseCountSensVsChrMap
#SBATCH --mail-type=END,FAIL
#SBATCH --time=4:59:59
#SBATCH --constraint=""

# goal
# -----
# Given a pause file and an expected pause count per read, per base bed file, plot the pause count and expected pause
# count per genomic window across all contigs in a given fasta fai file.

# usage
#------
# sbatch run_plot_pause_count_sense_vs_chr_map.sh $pause_file $pause_prob_per_read_per_single_base_bed \
#   $fasta_fai $output_dir
# NOTE: \ means the command continues on the next line.
# $pause_file: pause file in the standard pause format used by our pipeline (get_sgm_dnascent_pauses.sh)
# $pause_prob_per_read_per_single_base_bed: bed file with the expected pause count per read, per base as made
#    by run_get_sgm_pause_sensitivities.sh
# $fasta_fai: fasta fai file of a reference genome (in fasta format)
# $output_dir: directory to save the output plots/data files

# outputs
# --------
# A few files are sent to the output directory:
# - a count of pauses per genomic window across all contigs in the fasta fai file.
# - a count of expected pauses per genomic window across all contigs in the fasta fai file.
# - a plot of the pause count and expected pause count per genomic window across all contigs in the fasta fai file.

# stop execution if any command fails
set -e

# load python
source load_package.sh -python -R -imagemagick

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# assign arguments to variables
pause_file=$1
pause_prob_per_read_per_single_base_bed=$2
fasta_fai=$3
output_dir=$4

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
  echo >&2 "ERROR: incorrect number of arguments provided"
  echo >&2 "USAGE: sbatch <program_name.sh> pause_file pause_prob_per_read_per_single_base_bed fasta_fai output_dir"
  echo >&2 "pause_file: pause file in the standard pause format used by our pipeline (get_sgm_dnascent_pauses.sh)"
  echo >&2 "pause_prob_per_read_per_single_base_bed: bed file with the expected pause count per read, per base as "
  echo >&2 "                                         made by run_get_sgm_pause_sensitivities.sh"
  echo >&2 "fasta_fai: fasta fai file of a reference genome (in fasta format)"
  echo >&2 "output_dir: directory to save the output plots/data files"
  echo >&2 "For more details, please see the documentation in the script. Exiting..."
  exit 1
fi

# check that the input files exist
if [ ! -f "$pause_file" ]; then
  echo >&2 "ERROR: pause file does not exist. Exiting..."
  exit 1
fi
if [ ! -f "$pause_prob_per_read_per_single_base_bed" ]; then
  echo >&2 "ERROR: pause_prob_per_read_per_single_base_bed does not exist. Exiting..."
  exit 1
fi
if [ ! -f "$fasta_fai" ]; then
  echo >&2 "ERROR: fasta fai file does not exist. Exiting..."
  exit 1
fi

# make the output directory if it does not exist
mkdir -p "$output_dir"

# set output files
pause_bedgraph="$output_dir/pause_count_per_genomic_window.bedgraph"
sens_bedgraph="$output_dir/expected_pause_count_per_genomic_window.bedgraph"

# create temporary file
tmp_bed=$(mktemp "$tmpDir"/tmp.XXXXXXXXXX.bed)

# get pause count per genomic window
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses > "$tmp_bed"
bash convert_bed_to_coverage.sh "$fasta_fai" 1000 "$tmp_bed" --invert > "$pause_bedgraph"

# get expected pause counts per genomic window
grep -v '^#' "$pause_prob_per_read_per_single_base_bed" | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $7, $6}' |\
    bash sum_per_read_per_base_bed_to_genome_window_bed.sh - "$fasta_fai" 1000 chrM > "$sens_bedgraph"

# plot counts and expected counts vs genomic coordinates
cd plotting_and_short_analyses || exit
Rscript plot_pause_count_sens_vs_chr_map.R "$fasta_fai" "$pause_bedgraph" "$sens_bedgraph" \
  "$output_dir"/pause_count_sens_vs_chr_map.svg

# find the svg file and convert it to png
convert -density 1000 "$output_dir"/pause_count_sens_vs_chr_map.svg \
  "$output_dir"/pause_count_sens_vs_chr_map.png

# remove temporary files
rm -rf "$tmpDir"