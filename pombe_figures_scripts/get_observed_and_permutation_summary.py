import sys
import glob
import os
import numpy as np
import pandas as pd
from collections import defaultdict

def count_distances_in_windows(distances):
    """Count distances in 1kb bins."""
    if len(distances) == 0:
        return {}
    
    max_dist = max(distances)
    bins = np.arange(0, max_dist + 1000, 1000)
    digitized = np.digitize(distances, bins) - 1  # -1 to make it 0-based
    
    # Create bin labels like "(0,1000]"
    bin_labels = [f"({bins[i]},{bins[i+1]}]" for i in range(len(bins)-1)]
    
    # Count occurrences in each bin
    counts = defaultdict(int)
    for d in digitized:
        if d < len(bin_labels):  # Ensure we don't go out of bounds
            counts[bin_labels[d]] += 1
    
    return counts

def process_permutation_files(permutation_pattern, num_permutations):
    """Process all permutation files and return counts DataFrame."""
    all_counts = []
    permutation_files = sorted(glob.glob(permutation_pattern))[:num_permutations]
    
    for i, file_path in enumerate(permutation_files, 1):
        try:
            # Read file and calculate distances
            data = pd.read_csv(file_path, sep='\t', header=None)
            distances = np.abs(data[1] - data[4])  
            
            # Count distances in bins
            counts = count_distances_in_windows(distances)
            all_counts.append(counts)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    # Get all possible bin names across all datasets
    all_bins = set()
    for counts in all_counts:
        all_bins.update(counts.keys())
    all_bins = sorted(all_bins, key=lambda x: float(x.split(',')[0][1:]))
    
    # Create filled DataFrame
    data = []
    for i, counts in enumerate(all_counts, 1):
        row = {'Dataset': f'Dataset_{i}'}
        for bin_name in all_bins:
            row[f'bin_{bin_name}'] = counts.get(bin_name, 0)
        data.append(row)
    
    df = pd.DataFrame(data)
    return df, all_bins

def calculate_summary_stats(counts_df, bin_columns):
    """Calculate mean and standard deviation for each bin."""
    summary = pd.DataFrame({
        'bin': [col.replace('bin_', '') for col in bin_columns],
        'permutation_mean': counts_df[bin_columns].mean().values,
        'permutation_sd': counts_df[bin_columns].std().values
    })
    return summary

def process_observed_file(observed_file, all_bins):
    """Process the observed distances file."""
    data = pd.read_csv(observed_file, sep='\t', header=None)
    distances = np.abs(data[1] - data[4])  
    counts = count_distances_in_windows(distances)
    
    # Create DataFrame with all bins (even if count is 0)
    obs_counts = []
    for bin_name in all_bins:
        obs_counts.append({'bin': bin_name, 'obs_count': counts.get(bin_name, 0)})
    
    return pd.DataFrame(obs_counts)

def main():
    if len(sys.argv) != 5:
        print("Usage: python get_observed_and_permutation_summary.py <observed_distances.txt> <permutation_pattern> <num_permutations> <output_dir>")
        sys.exit(1)
    
    observed_file = sys.argv[1]
    permutation_pattern = sys.argv[2]
    num_permutations = int(sys.argv[3])
    output_dir = sys.argv[4]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Process permutation files
    print("Processing permutation files...")
    counts_df, all_bins = process_permutation_files(permutation_pattern, num_permutations)
    bin_columns = [f'bin_{b}' for b in all_bins]
    
    # Save permutation counts
    counts_df.to_csv(os.path.join(output_dir, 'permutation_counts.csv'), index=False)
    
    # Calculate summary stats
    summary_df = calculate_summary_stats(counts_df, bin_columns)
    
    # 2. Process observed file
    print("Processing observed distances...")
    obs_df = process_observed_file(observed_file, all_bins)
    
    # 3. Merge observed and null stats
    merged_df = pd.merge(obs_df, summary_df, on='bin', how='left')
    
    # Save results
    merged_df.to_csv(os.path.join(output_dir, 'observed_vs_permutation_summary.csv'), index=False)
    print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
