#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J count_origins
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# Goal: To count number of origins (given a list of mid-point origin locations) per sliding window.

# Usage
# sbatch run_count_origins_per_window.sh <genome_fai_file> <origins_file> <window_size> <slide_size> <output_file>
## genome_fai_file: fasta.fai file with the contig and coordinates
## origins_file: text file with contig in the first column and origin midpoint in the second column.
## window_size: size of window (e.g. 1000 nucleotides)
## sliding_size: size of sliding window step (e.g. 1 nucleotide)
## output_file: text file with tab-separated columns: contig, start_window, end_window, and origin_count

# Load python
source load_package.sh -python

# Print the arguments for debugging purposes
echo "FASTA index file: $1"
echo "Origins file: $2"
echo "Window size: $3"
echo "Slide size: $4"
echo "Output file: $5"

# Run the Python script with user-specified arguments
python count_origins_per_window.py $1 $2 $3 $4 $5


