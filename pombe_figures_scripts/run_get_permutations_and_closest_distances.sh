#!/bin/bash

#SBATCH --mem-per-cpu=1G
#SBATCH -c 2
#SBATCH -p ei-short
#SBATCH -J permutation
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1:29:59


# Goal: to permutate initiation site intervals from an origin_forkSense.bed file, and to calculate the closest distance between the midpoint
# of the initiation site intervals from each permutation and the initiation site midpoint from a file of interest (e.g., Pu-seq origins).

# Usage: run_get_permutations_and_closest_distances.sh <origins_forkSense.bed> <reference_genome.fai> <num_permutations> <origins_other_study.txt> <output_directory>

# Stop execution if any command fails
set -e

# Load required packages
source load_package.sh -bedtools -python

# Load configuration
source config.sh

# Ensure the user provides five arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <origins_forkSense.bed> <reference_genome.fasta.fai> <num_permutations> <origins_other_study.txt> <output_directory>"
    exit 1
fi

# Assign the arguments to variables
origins_forkSense="$1"
reference_genome="$2"
num_permutations="$3"
origins_other_study="$4"
output_dir="$5"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$output_dir"/distances_permutations

# Process the origins file
echo "Processing origins file..."
awk -F' ' 'NR > 1 {print $1 "\t" $2 "\t" $3}' "$origins_forkSense" | sort -k1,1 -k2,2n > "$output_dir/temp_3col.bed"

# Map chromosome names and sort
awk 'BEGIN {
    map["NC_003424.3"] = "chr1"
    map["NC_003423.3"] = "chr2"
    map["NC_003421.2"] = "chr3"
    map["NC_001326.1"] = "chrM"
    FS = OFS = "\t"
}
{
    if ($1 in map) { $1 = map[$1] }
    print
}' "$output_dir/temp_3col.bed" | sort -k1,1 -k2,2n > "$output_dir/temp_3col_chrom.bed"

# Calculate midpoints and remove chrM
awk -F'\t' -v OFS='\t' '{
    if ($1 == "chrM") next
    midpoint = int(($2 + $3) / 2)
    print $1, midpoint, midpoint
}' "$output_dir/temp_3col_chrom.bed" | sort -k1,1 -k2,2n > "$output_dir/processed_origins_forkSense.bed"

# Get observed distances
echo "Calculating observed distances..."
bedtools closest -a "$output_dir/processed_origins_forkSense.bed" -b "$origins_other_study" > "$output_dir/observed_distances.txt"

# Run permutations
echo "Running $num_permutations permutations..."
bash get_permutations.sh "$output_dir/processed_origins_forkSense.bed" "$reference_genome" "$num_permutations" "$output_dir"

# Calculate closest distance between permutations and the file of interest
for ((i=1; i<=$num_permutations; i++)); do
    # Run bedtools closest for each permutation against the reference file
    output_distance_file="${output_dir}/distances_permutations/permutation_midpoint${i}_distance.txt"
    echo "Running bedtools closest for permutation $i..."
    bedtools closest -a "${output_dir}/permutations/permutation_midpoint${i}.bed" -b "$origins_other_study" > "$output_distance_file"
    echo "Finished bedtools closest for permutation $i."
done

python get_observed_and_permutation_summary.py "$output_dir/observed_distances.txt" "$output_dir/distances_permutations/*.txt" "$num_permutations" "$output_dir"


# Clean up temporary files
rm "$output_dir/temp_3col.bed" "$output_dir/temp_3col_chrom.bed"

echo "All done!"
