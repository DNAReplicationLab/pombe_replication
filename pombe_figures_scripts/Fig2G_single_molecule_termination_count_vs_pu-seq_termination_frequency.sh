#!/bin/bash

# Goal: Count single-molecule termination midpoints in 1 kb sliding windows (with 300 bp slide),
# then compare counts to Pu-seq termination frequency and plot the relationship.

# Usage: bash Fig2G_single_molecule_termination_count_vs_pu-seq_termination_frequency.sh /path/to/input_dir /path/to/output_dir

# Before using the code, please download the "pu-seq_termination_dataset.bed" (accession GSE62108) and the "pombe_reference_genome.fna" (ASM294v2). 
# Save those files in /path/to/input_dir

set -euo pipefail  

# Input and output directories
input_dir="$1"
output_dir="$2"

# Create output subdirectory
mkdir -p "$output_dir/Fig2G"

# Load required tools
source load_package.sh -R -samtools

# Index the reference genome for window generation
samtools faidx "$input_dir/pombe_reference_genome.fna"

# Step 1: Compute termination midpoint for each region in the BED file
awk -v OFS='\t' '{print $1, int(($2+$3)/2 + 0.9999)}' "$input_dir/terminations-IDS_65-R9-ASM294v2.bed" > "$output_dir/Fig2G/terminations_midpoint.txt"

# Step 2: Convert contig names to chr1/chr2/chr3 and remove mitochondrial chromosome
input_file="$output_dir/Fig2G/terminations_midpoint.txt"
output_file="$output_dir/Fig2G/terminations_midpoint_chr.txt"

awk 'BEGIN {
    # Mapping from NCBI contig names to standard chr names
    map["NC_003424.3"] = "chr1"
    map["NC_003423.3"] = "chr2"
    map["NC_003421.2"] = "chr3"
    FS = "\t"; OFS = "\t"
}
# Exclude mitochondrial chromosome
$1 == "NC_001326.1" { next }
# Convert contig name if it exists in the map
{
    if ($1 in map) {
        $1 = map[$1]
    }
    print
}' "$input_file" | sort -k1,1 -k2,2n > "$output_file"

# Step 3: Count terminations in 1000-bp windows with 300-bp slide
bash run_count_origins_per_window.sh \
  "$input_dir/pombe_reference_genome.fna.fai" \
  "$output_file" \
  1000 \
  300 \
  "$output_dir/Fig2G/count_terminations_per_sliding_window_1000ntwindow_300ntslide.txt"

# Step 4: Run R script to plot termination count vs Pu-seq termination frequency
Rscript plot_termination_count_vs_frequency.R \
  "$output_dir/Fig2G/count_terminations_per_sliding_window_1000ntwindow_300ntslide.txt" \
  "$input_dir/pu-seq_termination_dataset.bed" \
  "$output_dir"
