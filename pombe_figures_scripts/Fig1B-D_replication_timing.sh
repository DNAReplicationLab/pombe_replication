#!/bin/bash

# Goal: To generate replication timing profiles for Fig 1B–1D using plot_replication_timing.R
# Usage: bash Fig1B-D_replication_timing.sh /path/to/input_dir /path/to/output_dir

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
bed_cerevisiae="$input_dir/replication_time_s_cerevisiae.bed"

# Validate input files
for bed_file in "$bed_pombe" "$bed_cerevisiae"; do
  if [ ! -f "$bed_file" ]; then
    echo "Error: Input BED file not found: $bed_file"
    exit 1
  fi
done

# Ensure output directory exists
mkdir -p "$output_dir"

# Output plots
out_fig1b="$output_dir/Fig1B_replication_timing.pdf"
out_fig1c="$output_dir/Fig1C_replication_timing.pdf"
out_fig1d="$output_dir/Fig1D_replication_timing.pdf"

# Source R
source load_package.sh -R

# Generate plots
echo "Generating Figure 1B (S. pombe chrIII:1637464–1732232)..."
Rscript "$script_dir/plot_replication_timing.R" \
  "$bed_pombe" "$out_fig1b" "NC_003421.2" "1637464:1732232"

echo "Generating Figure 1C (S. cerevisiae chrVIII:215468–349289)..."
Rscript "$script_dir/plot_replication_timing.R" \
  "$bed_cerevisiae" "$out_fig1c" "chrVIII" "215468:349289"

echo "Generating Figure 1D (S. pombe chrIII:1641285–1715967)..."
Rscript "$script_dir/plot_replication_timing.R" \
  "$bed_pombe" "$out_fig1d" "NC_003421.2" "1641285:1715967"

echo "All plots generated in: $output_dir"
