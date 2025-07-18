#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Fig 4B using plot_rainplot.sh
# Usage: bash Fig4B_rainplots.sh /path/to/input_dir /path/to/output_dir 

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
bam_ids_65="$input_dir/IDS_65-R9-AW2224_assembly.mod.bam"
read_id_IDS_65_a="5081ecbb-dba7-4209-af97-a773782e2120"
read_id_IDS_65_b="fe381ba4-e800-4441-9b02-0d84fd5c1cc6"
read_id_IDS_65_c="e6ce34c5-9d50-414b-8a8b-cdfbf1d72d14"
read_id_IDS_65_d="6e1cb823-56ef-4f23-ad46-b1295cd193d4"

# Validate input files
for bam_file in "$bam_ids_65"; do
  if [ ! -f "$bam_file" ]; then
    echo "Error: Input BAM file not found: $bam_file"
    exit 1
  fi
done

# Parameters
thymidine_window="300"

# Output directories
out_fig4b_a="$output_dir/Fig4B_a_rainplot"
out_fig4b_b="$output_dir/Fig4B_b_rainplot"
out_fig4b_c="$output_dir/Fig4B_c_rainplot"
out_fig4b_d="$output_dir/Fig4B_d_rainplot"


# Ensure output directories exist
mkdir -p "$out_fig4b_a" "$out_fig4b_b" "$out_fig4b_c" "$out_fig4b_d" 

# Generate rainplots
echo "Generating rainplot for Figure 4B a..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_a" "$thymidine_window" "$out_fig4b_a"
echo "Figure 4B a complete."

echo "Generating rainplot for Figure 4B b..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_b" "$thymidine_window" "$out_fig4b_b"
echo "Figure 4B b complete."

echo "Generating rainplot for Figure 4B c..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_c" "$thymidine_window" "$out_fig4b_c"
echo "Figure 4B c complete."

echo "Generating rainplot for Figure 4B d..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_d" "$thymidine_window" "$out_fig4b_d"
echo "Figure 4B d complete."

