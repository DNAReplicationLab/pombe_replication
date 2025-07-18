#!/bin/bash

#SBATCH --mem=5G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotFltNasReads
#SBATCH --mail-type=END,FAIL
#SBATCH --time 4:59:59

# goal of program
# =================
# Select a few nascent reads at random that satisfy some criteria and produce one plot and one data file per read.
# Each plot shows raw and windowed analogue probabilities versus reference coordinate.
# Prerequisite: Run main pipeline of fast5 -> dnascent before this.

# program sample usage and input parameter description
# ====================================================
# sbatch <program.sh> n forkLenThres alignLenThres mod_bam mod_bam_left mod_bam_right forkSense_dir output_dir
# n: number of reads requested. WARNING: n should be small, like 10 - 100.
#    Do not use larger numbers unless you know what you are doing.
# forkLenThres: (in kb) only consider reads where forks of at least this length are detected.
# alignLenThres: (in kb) only consider reads with alignments of at least this length.
# mod_bam: BAM file containing modification data. will have been produced by pipeline
# mod_bam_left: BAM file containing probabilities of a left fork per thymidine position per read.
#   will have been produced by pipeline
# mod_bam_right: BAM file containing probabilities of a right fork per thymidine position per read.
#   will have been produced by pipeline
# forkSense_dir: directory containing fork sense data. will have been produced by pipeline.
# op_dir: directory where the n plots go

# stop execution if any command fails
set -e

# load variables from the command line
n_reads=$1;
fork_length_threshold=${2:-0};
alignment_length_threshold=$3
mod_bam=$4;
mod_bam_left=$5;
mod_bam_right=$6;
forksense_dir=$7;
op_dir=$8;

# check if eight arguments were passed
if [ $# -ne 8 ]; then
    echo "ERROR: Eight arguments are required. You passed $# arguments.";
    exit 1;
fi

# check if op_dir exists, if it does then exit the program
if [ -d "$op_dir" ]; then
    echo "Output directory $op_dir already exists. Exiting program.";
    exit 1;
else
    mkdir -p "$op_dir";
    op_dir=$(readlink -f "$op_dir");
fi

# check left and right fork files exist
left_fork_file="$forksense_dir"/leftForks_DNAscent_forkSense.bed
right_fork_file="$forksense_dir"/rightForks_DNAscent_forkSense.bed

if [[ ! -f "$left_fork_file" ]] || [[ ! -f "$right_fork_file" ]] || [[ ! -f "$mod_bam_left" ]] ||\
  [[ ! -f "$mod_bam_right" ]] || [ "$fork_length_threshold" -le 0 ]; then
    echo "fork information is not set properly, so proceeding without fork information";
    fork_set=0;
else
    fork_set=1;
fi

# multiply fork_length_threshold and alignment_length_threshold by 1000 to convert kb to b
fork_length_threshold=$((fork_length_threshold * 1000))
alignment_length_threshold=$((alignment_length_threshold * 1000))

# make a temporary file to store the read ids
read_id_file=$(mktemp -p "$op_dir");

# if fork information is set, then use it to filter the reads.
# concatenate left and right fork files, and filter lines such that the
# difference of third and second column must be at least fork_length_threshold, and the difference of the
# seventh and sixth columns must be at least alignment_length_threshold, and the first column must not be
# equal to chrM. Then, select n_reads unique read ids (value of fourth column) at random and send it to read_id_file.
if [ "$fork_set" -eq 1 ]; then
    cat "$left_fork_file" "$right_fork_file" |\
    awk -v fork_length_threshold="$fork_length_threshold" -v alignment_length_threshold="$alignment_length_threshold" \
    '{if ($3-$2 >= fork_length_threshold && $7-$6 >= alignment_length_threshold && $1 != "chrM") print $4}' |\
    sort | uniq | shuf -n "$n_reads" > "$read_id_file";
else
    # if fork information is not set, then just select n_reads unique read ids at random from the bam file
    # where the read length is at least alignment_length_threshold.
    samtools view -e "rlen>=$alignment_length_threshold"  "$mod_bam" | cut -f 1 | sort | uniq |\
      shuf -n "$n_reads" > "$read_id_file";
fi

# make the plots and collate into a pdf
sbatch --job-name=lpr plot_n_reads.sh "$read_id_file" "$mod_bam" "$mod_bam_left" "$mod_bam_right" "$forksense_dir"\
  /dev/null /dev/null "$op_dir"

# delete read id file after plotting has finished.
# we use the wait flag so that the parent launcher job sticks around till all plots have been made.
# the wait may not be necessary depending on the context, but it's nice to have.
sbatch --dependency=singleton --job-name=lpr --wait -p ei-short --wrap="rm $read_id_file;"