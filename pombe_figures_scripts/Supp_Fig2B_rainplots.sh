#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Supp Fig 2B using plot_rainplot.sh
# Usage: bash Fig1B-D_rainplots.sh /path/to/input_dir /path/to/output_dir 

set -euo pipefail

# Usage message
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_directory> <output_directory>"
  exit 1
fi

# Set input and output directories from user arguments
input_dir="$1"
output_dir="$2"

# Set script directory
script_dir="$(dirname "$0")"
plot_script="$script_dir/plot_rainplot.sh"

# Input BAM files and read IDs
bam_ids_47_r10="$input_dir/IDS_47-R10-ASM294v2.mod.bam"
read_id_ids_47_r10="a5eca2e9-ba16-44c2-88a1-4b4aa59abd2a"

bam_ids_65_r10="$input_dir/IDS_65-R10-ASM294v2.mod.bam"
read_id_ids_65_r10="354e890b-58b0-488d-8c65-e576a544e401"

# Validate input files
for bam_file in "$bam_ids_47_r10" "$bam_ids_65_r10"; do
  if [ ! -f "$bam_file" ]; then
    echo "Error: Input BAM file not found: $bam_file"
    exit 1
  fi
done

# Parameters
thymidine_window="300"

# Output directories
out_supp_fig2b_left="$output_dir/Supp_Fig2B_left_rainplot"
out_supp_fig2b_right="$output_dir/Supp_Fig2B_right_rainplot"


# Ensure output directories exist
mkdir -p "$out_supp_fig2b_left" "$out_supp_fig2b_right"

# Generate rainplots
echo "Generating rainplot for Supp Figure 2B left..."
bash "$plot_script" "$bam_ids_47_r10" "$read_id_ids_47_r10" "$thymidine_window" "$out_supp_fig2b_left"
echo "Supp Figure 2B left complete."

echo "Generating rainplot for Supp Figure 2B right..."
bash "$plot_script" "$bam_ids_65_r10" "$read_id_ids_65_r10" "$thymidine_window" "$out_supp_fig2b_right"
echo "Supp Figure 2B right complete."
