#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --ntasks 15
#SBATCH -p ei-medium
#SBATCH -J makeSubset
#SBATCH --mail-type=END,FAIL
#SBATCH --time=6:00:00

# load packages
source load_package.sh -samtools

# stop execution if any step fails
set -e

# inputs are: set of read ids requested, and old and new filenames
#   read ids should be in one column with no header
#   NOTE that a bam file could also mean a mod bam file.
read_list=
alignment_old=
alignment_new=
seq_sum_old=
seq_sum_new=
fork_sense_overall_old_dir=
fork_sense_overall_new_dir=

# extract a subset of bam file
if [ -f "$alignment_old" ]; then
  samtools view -@ 15 -b -N "$read_list" -h "$alignment_old" > "$alignment_new"_temp
fi

# extract a subset of the sequencing summary file
if [ -f "$seq_sum_old" ]; then
  {
    head -n 1 "$seq_sum_old";
    grep -F -f "$read_list" "$seq_sum_old";
  } > "$seq_sum_new"
fi

if [ -d "$fork_sense_overall_old_dir" ]; then
  # make the fork sense directory if need be
  mkdir -p "$fork_sense_overall_new_dir"

  prefixes=('rightForks' 'leftForks' 'origins' 'terminations')

  for prefix in "${prefixes[@]}"
  do
    filename="${prefix}"_DNAscent_forkSense.bed
    grep -F -f "$read_list" "$fork_sense_overall_old_dir"/"$filename"  > "$fork_sense_overall_new_dir"/"$filename"
  done
fi

# sort and index the new bam files if they exist
if [ -f "$alignment_new"_temp ]; then
  samtools sort -@ 15 -o "$alignment_new" "$alignment_new"_temp
  samtools index -@ 15 "$alignment_new"
  rm "$alignment_new"_temp;
fi