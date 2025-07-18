#!/bin/bash

# Goal: Count initiation midpoints within 1 kb from Pu-seq origin, group by efficiency, and plot.
# Usage: bash Fig2F_single_molecule_initiation_count_vs_pu-seq_efficiency.sh /path/to/input_dir /path/to/output_dir

# Before using the code, please download the "pu-seq_origin_dataset.bed" (accession GSE62108). 
# Save the file in /path/to/input_dir

set -euo pipefail

# Get user-defined input/output directories
input_dir="$1"
output_dir="$2"

# Create output subdirectory
mkdir -p "$output_dir/Fig2F"

# Load required tools
source load_package.sh -R

# Run R script
Rscript plot_count_vs_efficiency.R "$input_dir/initiations-IDS_65-R9-ASM294v2.bed" "$input_dir/pu-seq_origin_dataset.bed" "$output_dir"
