#!/bin/bash

# Goal: To generate replication timing profiles for Supp Fig 2B using plot_replication_timing.R
# Usage: bash Supp_Fig2B_replication_timing.sh /path/to/input_dir /path/to/output_dir

set -euo pipefail

# Usage message
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_directory> <output_directory>"
  exit 1
fi

# User-specified directories
input_dir="$1"
output_dir="$2"

# Script directory
script_dir="$(dirname "$0")"

# Input BED files
bed_pombe="$input_dir/replication_time_s_pombe.bed"

# Validate input files
for bed_file in "$bed_pombe"; do
  if [ ! -f "$bed_file" ]; then
    echo "Error: Input BED file not found: $bed_file"
    exit 1
  fi
done

# Ensure output directory exists
mkdir -p "$output_dir"

# Output plots
out_supp_fig2b_left="$output_dir/Supp_Fig2B_left_replication_timing.pdf"
out_supp_fig2b_right="$output_dir/Supp_Fig2B_right_replication_timing.pdf"

# Source R
source load_package.sh -R

# Generate plots
echo "Generating Supp Figure 2B left (S. pombe NC_003423.3:1662247–1746964)..."
Rscript "$script_dir/plot_replication_timing.R" \
  "$bed_pombe" "$out_supp_fig2b_left" "NC_003423.3" "1662247:1746964"

echo "Generating Supp Figure 2B right (S. pombe NC_003421.2:1893652–1994859)..."
Rscript "$script_dir/plot_replication_timing.R" \
  "$bed_pombe" "$out_supp_fig2b_right" "NC_003421.2" "1893652:1994859"

echo "All plots generated in: $output_dir"
