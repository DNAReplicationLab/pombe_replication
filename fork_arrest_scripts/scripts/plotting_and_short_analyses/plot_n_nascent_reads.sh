#!/bin/bash

#SBATCH --mem=5G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotNasReads
#SBATCH --mail-type=END,FAIL
#SBATCH --time 4:59:59

# goal of program
# =================
# Select a few nascent reads at random and produce one plot and one data file per read.
# Each plot shows raw and windowed analogue probabilities versus reference coordinate.
# Prerequisite: Run main pipeline of fast5 -> dnascent before this.

# program sample usage and input parameter description
# ====================================================
# sbatch <program.sh> n mod_bam_file brdu_data_file op_dir win_size threshold [mod_threshold]
# n: number of reads requested. WARNING: n should be small, like 10 - 100.
#    Do not use larger numbers unless you know what you are doing.
# mod_bam_file: BAM file containing modification data. will have been produced by pipeline
# brdu_data_file: tab separated without headers, first column read id, second column mean brdu density across read.
#                 Have a look in the pipeline output folder for the file.
# op_dir: directory where the n plots go, we recommend an empty directory unless you know what you are doing.
# win_size: window size (in number of thymidines) for windowing analogue probabilities
# thres: threshold above which reads are marked as nascent, used to select reads from brdu_data_file
# mod_threshold: (optional, default 0.5) threshold above which a thymidine is called as BrdU in the windowed data.
#               [] means the argument is optional. If you want to specify the argument, remove the brackets and
#               write the argument. If you don't want to specify the argument, just don't write anything.

# stop execution if any command fails
set -e

# load variables from the command line
n_reads=$1;
mod_bam=$2;
brdu_data_file=$3;
op_dir=$4;
win_size=$5;
thres=$6;
mod_threshold=${7:-0.5};

# make the output directory if needed
mkdir -p "$op_dir"

# select reads at random
< "$brdu_data_file" awk '{if($2 > '"$thres"'){print $1}}' | shuf | head -n "$n_reads" > "$op_dir"/read_ids.txt

# check that all reads are unique; it is possible for read ids to repeat in the input file.
# We cannot deal with this situation right now, so we print an error and exit.
n_repeated_reads=$(sort "${op_dir}/read_ids.txt" | uniq -d | wc -l)
if [ "$n_repeated_reads" -ne 0 ]
then
  echo "Read ids in the input file are not unique. Exiting."
  exit 1;
fi

# load samtools
current_dir=$(pwd)
main_dir=$(cd .. && pwd)
cd "$main_dir"
source load_package.sh -samtools
cd "$current_dir"

# subset modbam file to only include these reads
mod_bam_subset="${op_dir}/mod_bam_subset.bam"
samtools view -b -N "$op_dir"/read_ids.txt "$mod_bam" > "${mod_bam_subset}.tmp"
samtools sort -o "$mod_bam_subset" "${mod_bam_subset}.tmp"
samtools index "$mod_bam_subset"
rm "${mod_bam_subset}.tmp"

# plot the selected reads
# first, we send the plot jobs to a temporary file
submitted_jobs="${op_dir}/submitted_jobs"
< "$op_dir"/read_ids.txt \
    awk '{print "bash limited_plot_read.sh -t '"$mod_threshold"' '"$mod_bam_subset"' " $0 " '"$win_size"' '"$op_dir"' "}' \
    > "$submitted_jobs"

cd "$main_dir"
launchJob=$(sbatch -c 1 --mem=10G --time=1:00:00 -p ei-medium --array=1-"$n_reads" -J lPlotReadArray --mail-user= \
  -D "$current_dir" job_array_converter.sh "$submitted_jobs")
jid=${launchJob##* }
cd "$current_dir"

# launch a dummy job so that we wait till plot jobs have finished
sbatch --dependency=afterany:"$jid" --job-name=sleep_10 --wait -p ei-short --wrap="sleep 10;"

# print status of plotting jobs
sacct -j "$jid" --format=JobID,State,ExitCode,MaxRSS,Elapsed,Start,End