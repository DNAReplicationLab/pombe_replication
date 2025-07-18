#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 6
#SBATCH -p ei-medium
#SBATCH -J heatmap
#SBATCH --mail-type=START,END,FAIL
#SBATCH --time 23:59:59

# Goal: to extract raw and windowed data from reads within a specified genomic location
# and to plot it.

# Usage: sbatch modbam_to_heatmap.sh <mod.bam> <contig> <start> <end> <mod_code>
# <position> <threshold> <window_size> <output_dir>

# Input file: mod.bam

# Genomic location:
#  contig: the contig (e.g., chromosome) of the region of interest
#  start: the start position of the region of interest
#  end: the end position of the region of interest

# Modification code: mod_code
#  DNA base modification of interest
#  e.g., T for any T modification, h for 5hmC, m for 5mC, a for 6mA

# Position of DNA base modification: position
#  must be either not_use_ref or use_ref.
#  not_use_ref: the position is the value of the column forward_read_position
#  use_ref: the position is the value of the column ref_position

# Threshold for modification calling: threshold
#  e.g., a thres = 0.5 means that a base with a modification probability >= 0.5 is modified
#  and < 0.5 is unmodified

# Window size: window_size
#  number of bases that are used to window the data after thresholding
#  e.g., a win = 300 means that data is windowed in 300 bases after thresholding

# Output directory: output_dir

# Output files:
#  heatmap_$contig_$start-$end.svg: a heatmap with the DNA base modification "mod_qual" data.
#  heatmap_$contig_$start-$end.tsv: a csv file with the modification information from
#  the reads in the specified region.
#  read_ids_${contig}_${start}-${end}_ordered.txt: a list with read_ids in the same order as 
#  in the heatmap.

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <input_bam> <contig> <start> <end> <mod_code> <position> <threshold> <window_size> <output_dir>"
    exit 1
fi

# Define variables
input_bam=$1
contig=$2
start=$3
end=$4
mod_code=$5
position=$6
threshold=$7
window_size=$8
output_dir=$9

# Create output directory
mkdir -p "$output_dir" || { echo "Failed to create output directory: $output_dir"; exit 1; }

# Create a temporary directory for storing intermediate files
temp_dir=$(mktemp -d -p "$output_dir" tempdir.XXXXXX) || { echo "Failed to create temporary directory"; exit 1; }
echo "Temporary directory created: $temp_dir"

# Function to check if a command succeeded
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: Command failed - $1"
        rm -rf "$temp_dir"
        exit 1
    fi
}

# Load necessary packages
source load_package.sh -R -python -samtools -modkit
check_command "Loading required packages"

# Calculate dynamic chunk size based on BAM file size
bam_size=$(stat -c%s "$input_bam")
chunk_size=$((bam_size / 1000)) # Adjust divisor to control chunk size
chunk_size=$((chunk_size < 10000 ? 10000 : chunk_size)) # Ensure chunk size is at least 10,000
chunk_size=$(( (chunk_size + 999) / 1000 * 1000 )) # Round up to nearest multiple of 1000
echo "Using chunk size: $chunk_size"

# Define the final output files
combined_tsv="$output_dir/heatmap_${contig}_${start}-${end}.tsv"
combined_svg="$output_dir/heatmap_${contig}_${start}-${end}.svg"
read_ids_file="$output_dir/read_ids_${contig}_${start}-${end}_ordered.txt"

# Initialize the combined TSV file
> "$combined_tsv"
> "$read_ids_file"

# Process the BAM file in chunks
current_start=$start
while [ "$current_start" -lt "$end" ]; do
    current_end=$((current_start + chunk_size - 1))
    [ "$current_end" -gt "$end" ] && current_end=$end

    chunk_output_dir="$temp_dir/chunk_${current_start}-${current_end}"
    mkdir -p "$chunk_output_dir" || { echo "Failed to create chunk directory: $chunk_output_dir"; exit 1; }

    region_bam="$chunk_output_dir/heatmap.bam"
    region_bam_ts="$chunk_output_dir/heatmap_alldata.tsv"
    raw_detect_ts="$chunk_output_dir/heatmap_rawdata.tsv"
    win_detect_ts="$chunk_output_dir/heatmap_windata.tsv"

    # Extract reads from the specified region
    samtools view -b "$input_bam" "$contig:$current_start-$current_end" > "$region_bam"
    check_command "samtools view for region $current_start-$current_end"

    samtools index "$region_bam"
    check_command "samtools index for region $current_start-$current_end"

    # Extract modification data using modkit
    modkit extract full "$region_bam" "$region_bam_ts"
    check_command "modkit extract for region $current_start-$current_end"

    # Extract raw modification data
    python extract_raw_mod_data.py "$mod_code" "$position" "$region_bam_ts" "$raw_detect_ts"
    check_command "extract_raw_mod_data.py for region $current_start-$current_end"

    # Window the modification data
    python window_mod_data.py "$threshold" "$window_size" "$current_start" "$current_end" "$raw_detect_ts" "$win_detect_ts"
    check_command "window_mod_data.py for region $current_start-$current_end"

    # Append windowed data to the combined TSV file
    [ -f "$win_detect_ts" ] && cat "$win_detect_ts" >> "$combined_tsv"
    check_command "Appending windowed data for region $current_start-$current_end"

    echo "Processed region $current_start-$current_end"
    current_start=$((current_end + 1))
done

# Generate the combined heatmap
Rscript ./plotting_and_short_analyses/plot_heatmap_pombe_manuscript.R "$combined_tsv" "$combined_svg" "$read_ids_file"
check_command "Rscript plot_heatmap_pombe_manuscript.R"

echo "Combined heatmap generated: $combined_svg"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$temp_dir"
check_command "Cleaning up temporary files"

echo "Processing complete."
