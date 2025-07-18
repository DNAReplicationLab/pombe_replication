#!/bin/bash

# Goal: To calculate the fraction of left forks in 1 kb windows and plot it against Pu-seq (Daigaku et al 2015) fraction of left forks.
# Usage: bash Fig2D_single_molecule_leftward_fork_vs_pu-seq.sh /path/to/input_dir /path/to/output_dir

# Before using the code, please download the "pu-seq_fork_dataset.bed" (accession GSE62108) and the "pombe_reference_genome.fna" (ASM294v2). 
# Save those files in /path/to/input_dir

set -euo pipefail

# Get user-defined input/output directories
input_dir="$1"
output_dir="$2"

# Create output subdirectory
mkdir -p "$output_dir/Fig2D"

# Load required tools
source load_package.sh -bedtools -samtools -R

# Index the reference genome
samtools faidx "$input_dir/pombe_reference_genome.fna"

# Extract first 3 columns from DNAscent fork files and sort
awk '{print $1, $2, $3}' OFS="\t" < "$input_dir/left_forks-IDS_65-R9-ASM294v2.bed" \
  | sort -k1,1 -k2,2n > "$output_dir/Fig2D/left_sorted.bed"

awk '{print $1, $2, $3}' OFS="\t" < "$input_dir/right_forks-IDS_65-R9-ASM294v2.bed" \
  | sort -k1,1 -k2,2n > "$output_dir/Fig2D/right_sorted.bed"

# Genome coverage
bedtools genomecov -i "$output_dir/Fig2D/left_sorted.bed" \
  -g "$input_dir/pombe_reference_genome.fna.fai" -d \
  > "$output_dir/Fig2D/left_genomecov.bed"

bedtools genomecov -i "$output_dir/Fig2D/right_sorted.bed" \
  -g "$input_dir/pombe_reference_genome.fna.fai" -d \
  > "$output_dir/Fig2D/right_genomecov.bed"

# Calculate averages per 1kb bin 
bash run_calculate_bins.sh "$output_dir/Fig2D/left_genomecov.bed" \
  "$output_dir/Fig2D/left_genomecov_average.bed"

bash run_calculate_bins.sh "$output_dir/Fig2D/right_genomecov.bed" \
  "$output_dir/Fig2D/right_genomecov_average.bed"

# Convert averages to bed with fork direction
left_file="$output_dir/Fig2D/left_genomecov_average.bed"
right_file="$output_dir/Fig2D/right_genomecov_average.bed"
left_final="$output_dir/Fig2D/left_average_final.bed"
right_final="$output_dir/Fig2D/right_average_final.bed"

awk 'NR > 1 {print $1, $2, $3, $5, "left"}' OFS="\t" "$left_file" > "$left_final"
awk 'NR > 1 {print $1, $2, $3, $5, "right"}' OFS="\t" "$right_file" > "$right_final"

echo -e "chromosome\tstart\tend\tcount\tfork_direction" | cat - "$left_final" > temp && mv temp "$left_final"
echo -e "chromosome\tstart\tend\tcount\tfork_direction" | cat - "$right_final" > temp && mv temp "$right_final"

# Merge and relabel chromosomes
merged_file="$output_dir/Fig2D/right_and_left_final.bed"
merged_chr_file="$output_dir/Fig2D/right_and_left_final_chr.bed"

cat "$left_final" "$right_final" > "$merged_file"

awk 'BEGIN {
    map["NC_003424.3"] = "chr1"
    map["NC_003423.3"] = "chr2"
    map["NC_003421.2"] = "chr3"
    map["NC_001326.1"] = "chrM"
    FS = OFS = "\t"
}
{
    if ($1 in map) $1 = map[$1]
    print
}' "$merged_file" | sort -k1,1 -k2,2n > "$merged_chr_file"

# Calculate fraction of leftward forks
fraction_file="$output_dir/Fig2D/fraction_leftward_forks.bed"

awk '
BEGIN {
    OFS="\t"
    print "chr", "start", "end", "left_fraction"
}
{
    key = $1 "\t" $2 "\t" $3
    if ($5 == "left") {
        left_counts[key] += $4
    } else if ($5 == "right") {
        right_counts[key] += $4
    }
}
END {
    for (key in left_counts) {
        right = right_counts[key] ? right_counts[key] : 0
        total = left_counts[key] + right
        fraction = left_counts[key] / total
        printf "%s\t%.2f\n", key, fraction
    }
}' "$merged_chr_file" | sort -k1,1 -k2,2n > "$fraction_file"

# Combine with Pu-seq data
puseq_file="$input_dir/pu-seq_fork_dataset.bed"
comparison_file="$output_dir/Fig2D/leftward_forks_fraction_dnascent_and_puseq.txt"

awk '
BEGIN {
    OFS="\t"
    print "chr", "start", "end", "left_fraction", "left_fraction_Puseq"
}
NR==FNR && FNR>1 {
    key = $1 "\t" $2 "\t" $3
    a[key] = $4
    keys[key] = 1
    next
}
FNR>1 && FNR != NR {
    key = $1 "\t" $2 "\t" $3
    b[key] = $4
}
END {
    for (key in keys) {
        split(key, parts, "\t")
        chr = parts[1]
        start = parts[2]
        end = parts[3]
        fracA = a[key]
        fracB = (key in b) ? b[key] : "N/A"
        print chr, start, end, fracA, fracB
    }
}
' "$fraction_file" "$puseq_file" | sort -k1,1 -k2,2n > "$comparison_file"

echo "Combined results saved to $comparison_file"

# Plot correlation
Rscript plot_correlation.R "$comparison_file" "$output_dir"

# Cleanup all intermediate files except for the final comparison and plot
rm -f \
  "$output_dir/Fig2D/left_sorted.bed" \
  "$output_dir/Fig2D/right_sorted.bed" \
  "$output_dir/Fig2D/left_genomecov.bed" \
  "$output_dir/Fig2D/right_genomecov.bed" \
  "$output_dir/Fig2D/left_genomecov_average.bed" \
  "$output_dir/Fig2D/right_genomecov_average.bed" \
  "$output_dir/Fig2D/left_average_final.bed" \
  "$output_dir/Fig2D/right_average_final.bed" \
  "$output_dir/Fig2D/right_and_left_final.bed" \
  "$output_dir/Fig2D/right_and_left_final_chr.bed" \
  "$output_dir/Fig2D/fraction_leftward_forks.bed"

echo "Cleanup complete. Final output saved:"
echo "  - Data: $comparison_file"
echo "  - Plot: $output_dir/Fig2D/Fig2D_leftward_fraction_DNAscent_vs_Puseq.svg"
