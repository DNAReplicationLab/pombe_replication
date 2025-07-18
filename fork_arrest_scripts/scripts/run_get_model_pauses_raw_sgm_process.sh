#!/bin/bash

#SBATCH --mem=4G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J getModelPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-50
#SBATCH --time 11:59:59
#SBATCH --constraint=""

# note: The calling script must mention how many jobs are needed in the job array
#       (job array means many jobs are run in parallel).
#       We've set the number of jobs as 51 (#SBATCH --array=0-50) as an example.
#       Each job is responsible for detecting pauses in a continuous N-line-block of the input file,
#       where N is of order 1000 (pls read the code to see what value is used for N).
#       So, if the input file has say 30,500 lines, then set the number of jobs as at least 31.
#       If more than 31 jobs are specified, then the extra jobs will fail and this is not a problem.

# goal
# ----
# use the cut-and-align procedure to detect pauses

# usage
# -----
# sbatch run_get_model_pauses_raw_sgm_process.sh ltForks rtForks modBam opFile low high width iter tol
# ltForks: left forks file, in forkSense format
# rtForks: right forks file, in forkSense format
# opFile: output file
# low: lowest level of sigmoid
# high: highest level of sigmoid
# width: width of sigmoid
# iter: (optional, default 20) maximum number of iterations used by numerical procedure
# tol: (optional, default 0.03) tolerance used by numerical procedure

# To understand what iter and tol mean, please consult the program that uses them as inputs.

# preamble
# --------

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load python
source load_package.sh -python -samtools

# get configuration
source config.sh

# check that at least seven inputs are present on the command line
if [ $# -lt 7 ]; then
  echo "usage: sbatch run_get_model_pauses_raw_sgm_process.sh ltForks rtForks modBam opFile low high width iter tol"
  echo "iter and tol are optional"
  exit 1
fi

# set directories and output file
ltForks=$1
rtForks=$2
bam=$3
opFile=$4
sigmoid_low=$5
sigmoid_high=$6
sigmoid_width=$7
numerical_iter=${8:-20}
numerical_tol=${9:-0.03}

# check that the input files exist
if [ ! -f "$ltForks" ] || [ ! -f "$rtForks" ] || [ ! -f "$bam" ]; then
  echo "input files do not exist"
  exit 1
fi

# check that the bam index file exists
if [ ! -f "$bam".bai ]; then
  echo "bam index file does not exist"
  exit 1
fi

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# generate output split across several files
pauseFile="$opFile"_part_"$SLURM_ARRAY_TASK_ID"

# set the script information
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# gather forks and split into 1000 forks per job and process
start_line=$(echo "$SLURM_ARRAY_TASK_ID"*1000 + 1 | bc)
end_line=$(echo "$SLURM_ARRAY_TASK_ID"*1000 + 1000 | bc)
{
    < "$ltForks" awk '{print $1 " " $2 " " $3 " " $4 " " $4 "_" $5 "_" $6 "_" $7 "_" $8 "_L_" $2 "_" $3}';
    < "$rtForks" awk '{print $1 " " $2 " " $3 " " $4 " " $4 "_" $5 "_" $6 "_" $7 "_" $8 "_R_" $2 "_" $3}';
} | sed -n "$start_line,$end_line p" |\
    sed '1i\contig start end read_id alt_read_id' > "$tmpDir"/fork_file

# extract read ids from the forks file
read_ids="$tmpDir"/read_ids
grep -v '^#' "$tmpDir"/fork_file | awk '{print $4}' | sort | uniq > "$read_ids"

# make a subset mod bam file with just the reads of interest
bam_subset="$tmpDir"/mod_bam_subset.bam
samtools view -b -N "$read_ids" "$bam" > "${bam_subset}.tmp"
samtools sort -o "$bam_subset" "${bam_subset}.tmp"
samtools index "$bam_subset"

{
  # output script information
  echo "$infoStr"

  < "$tmpDir"/fork_file \
      python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$bam_subset" |\
      sed '1i\detectIndex\tposOnRef\tprobBrdU' |\
      python get_model_pauses_raw_sgm.py --width "$sigmoid_width" --thres 0.5 \
         --low "$sigmoid_low" --high "$sigmoid_high" --iter "$numerical_iter" --tol "$numerical_tol" |\
      python process_model_pause_results_raw_sgm.py --width "$sigmoid_width" --low "$sigmoid_low" --high "$sigmoid_high"
} > "$pauseFile"

# remove temporary directory
rm -rf "$tmpDir"