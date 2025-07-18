#!/bin/bash

# Goal: To plot the distribution of distances
#       a) from observed DNAscent initiation midpoints to other datasets, and
#       b) from 1000 random permutations of the DNAscent initiation genomic intervals to other datasets.
# Usage: bash Supp_Fig3A-D_closest_distances.sh /path/to/input_dir /path/to/output_dir

# Before using the code, please download:
## the "origins_daigaku_et_al_2015.bed" (accession GSE62108) as a tab-separated file with three columns (chromosome, start, end)
## the "origins_segurado_et_al_2003.bed" (https://pombe.oridb.org/paper_data.php?id=14566325) as a tab-separated file with three columns (chromosome, start, end)
## the "origins_heichinger_et_al_2006.bed" (https://pombe.oridb.org/paper_data.php?id=17053780) as a tab-separated file with three columns (chromosome, start, end)
## the "origins_hayashi_et_al_2007.bed" (https://pombe.oridb.org/paper_data.php?id=17304213) as a tab-separated file with three columns (chromosome, start, end)
## the "pombe_reference_genome.fna" (ASM294v2).
# Save those files in /path/to/input_dir

set -euo pipefail

# Get user-defined input/output directories
input_dir="$1"
output_dir="$2"

# Create output subdirectory
mkdir -p "$output_dir/Supp_Fig3A-D"

# Load required tools
source load_package.sh -bedtools -samtools -R

# Index the reference genome
samtools faidx "$input_dir/pombe_reference_genome.fna"

# ---- Supp Fig 3A - Daigaku et al 2015 (Pu-seq) ----
bash run_get_permutations_and_closest_distances.sh \
  "$input_dir/initiations-IDS_65-R9-ASM294v2.bed" \
  "$input_dir/pombe_reference_genome.fna.fai" \
  1000 \
  "$input_dir/origins_daigaku_et_al_2015.bed" \
  "$output_dir/Supp_Fig3A-D/closest_distance_to_daigaku_et_al_2015"

# ---- Supp Fig 3B - Segurado et al 2003 ----
bash run_get_permutations_and_closest_distances.sh \
  "$input_dir/initiations-IDS_65-R9-ASM294v2.bed" \
  "$input_dir/pombe_reference_genome.fna.fai" \
  1000 \
  "$input_dir/origins_segurado_et_al_2003.bed" \
  "$output_dir/Supp_Fig3A-D/closest_distance_to_segurado_et_al_2003"

# ---- Supp Fig 3C - Heichinger et al 2006 ----
bash run_get_permutations_and_closest_distances.sh \
  "$input_dir/initiations-IDS_65-R9-ASM294v2.bed" \
  "$input_dir/pombe_reference_genome.fna.fai" \
  1000 \
  "$input_dir/origins_heichinger_et_al_2006.bed" \
  "$output_dir/Supp_Fig3A-D/closest_distance_to_heichinger_et_al_2006"

# ---- Supp Fig 3D - Hayashi et al 2007 ----
bash run_get_permutations_and_closest_distances.sh \
  "$input_dir/initiations-IDS_65-R9-ASM294v2.bed" \
  "$input_dir/pombe_reference_genome.fna.fai" \
  1000 \
  "$input_dir/origins_hayashi_et_al_2007.bed" \
  "$output_dir/Supp_Fig3A-D/closest_distance_to_hayashi_et_al_2007"

# Plot the distances as histograms
Rscript plot_histogram_distances_SuppFig3A-D.R "$output_dir/Supp_Fig3A-D"
