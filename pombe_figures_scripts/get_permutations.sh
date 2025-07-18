#!/bin/bash
#SBATCH --mem-per-cpu=2G
#SBATCH -c 2
#SBATCH -p ei-short
#SBATCH -J permutation
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1:59:59

# Usage: get_permutations.sh <input_bed> <genome_fai> <num_permutations> <output_dir>

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_bed> <genome_fai> <num_permutations> <output_dir>"
    exit 1
fi

source load_package.sh -bedtools

input_file="$1"
genome_file="$2"
num_permutations="$3"
output_dir="$4"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
mkdir -p "$output_dir"/permutations

# Loop for permutations
for ((i=1; i<=num_permutations; i++)); do
    # Run bedtools shuffle and name the output file as permutationX.bed
    output_file1="$output_dir/permutations/permutation${i}.bed"
    output_file2="$output_dir/permutations/permutation_midpoint${i}.bed"

    echo "Running permutation $i..."

    bedtools shuffle -i "$input_file" -g "$genome_file" > "$output_file1"

    # Transform the shuffled BED file
    awk '{midpoint = int(($2 + $3) / 2); print $1 "\t" midpoint "\t" midpoint}' "$output_file1" | \
    awk 'BEGIN {FS=OFS="\t"} {
        if ($1 == "NC_001326.1") next; 
        if ($1 ~ /^NC_003424.3/) { sub(/^NC_003424.3/, "chr1", $1) } 
        else if ($1 ~ /^NC_003423.3/) { sub(/^NC_003423.3/, "chr2", $1) } 
        else if ($1 ~ /^NC_003421.2/) { sub(/^NC_003421.2/, "chr3", $1) }
        print
    }' | \
    sort -k1,1 -k2,2n > "$output_file2"

    echo "Finished permutation $i."
done
