#!/bin/bash

#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J getModelFit
#SBATCH --mail-type=END,FAIL
#SBATCH --time 47:59:59

# goal
# -------------
# fit a model to (windowed) analogue probability data from many forks

# usage
# -----
# sbatch run_get_model_fit.sh $mod_bam $left_forks $right_forks $op_dir $tag
# can use bash in place of sbatch
# mod_bam: .mod.bam file with analogue modification probabilities
# left_forks: file with left forks in forksense format (see note). file must not have too many forks (see note).
# right_forks: file with right forks in forksense format (see note). file must not have too many forks (see note).
# op_dir: output directory
# tag: a short string to be used as a suffix in output file names

# note about forksense file format
# ---------------------------------
# space-separated with no header
# each row is: contig startFork endFork readID contig startAlignment endAlignment orn
# startFork always < endFork irrespective of whether its a left or right fork
# orn is fwd/rev
# this is the format of files output by DNAscent forkSense

# note about forksense file size
# ---------------------------------
# Script tested with each fork file with 100s of forks but not 1000s of forks.
# If you are going to use very large numbers of forks, then test the script to see if it can handle it,
# adjust memory etc.

# output files
# ------------
# two files are output: investigate_fit_params_<tag> and investigate_fit_<tag>
# the params file shows the best-fit parameters and the parameters within the 95% CI.
# NOTE: these parameters are only as good as the input-parameter-grid in get_model_fit.py
#       so adjust that if you want to scan over a coarser/finer grid, wider range etc.
# the investigate fit file shows coarse-grained data for each fork and how well it fits the model.
# a sample line and its meaning is shown below
# note that the sample line below is not real data and that the format may change over time if the code is changed.

# ce735f34-4705-430c-b38c-7d3962ae29ae_chrV_679828_745774_rev_L_681138_701477 0.64 700625 701471 36 0 0.609 NA NA good_short
# there are ten columns above.
# first column is forkName in the format 'readID_contig_startAlign_endAlign_orn_dirn_startFork_endFork'
# second column is mean BrdU
# third, fourth columns are the window boundaries over which the mean was recorded
# fifth column is the x value on the model at which that window of data fits best
# sixth column means whether the window is rejected (1 means true)
# seventh column is the goodness of fit of that window of data compared with model
# eight and ninth are the spread in the data, they are used only while reporting aggregates
# tenth column is a label that says if the fit is good or bad, if the fork is long or short etc.

# preamble
# --------

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load python, samtools
source load_package.sh -python -samtools

# load configuration variables
source config.sh

if [ "$#" -ne 5 ]
then
  echo "Incorrect number of arguments! "
  exit
fi

tag=$5
mkdir -p "$4";
fit_dir=$(cd "$4"; pwd)
mod_bam=$1
right_forks=$3
left_forks=$2

# set the script information
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# extract forks of interest
{
    < "$left_forks" awk 'BEGIN {OFS=" "}{print $2, $3, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_L_" $2 "_" $3, $4}'
    < "$right_forks" awk 'BEGIN {OFS=" "}{print $2, $3, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_R_" $2 "_" $3, $4}'
}  | sed '1i\start end alt_read_id read_id' | sed "1i$infoStr" > "$fit_dir"/forks_investigate_fit_"${tag}"

# extract read ids from the forks file
read_id_processed="$fit_dir"/read_id_processed_"${tag}"
grep -v '^#' "$fit_dir"/forks_investigate_fit_"${tag}" | awk '{print $4}' | sort | uniq > "$read_id_processed"

# make a subset mod bam file with just the reads of interest
bam_subset="$fit_dir"/mod_bam_subset_"${tag}".bam
samtools view -b -N "$read_id_processed" "$mod_bam" > "${bam_subset}.tmp"
samtools sort -o "$bam_subset" "${bam_subset}.tmp"
samtools index "$bam_subset"
rm "${bam_subset}.tmp"

# perform calculation
{
  echo "$infoStr"
  < "$fit_dir"/forks_investigate_fit_"${tag}" \
    python get_raw_data_from_modBAM.py --piped-regions --fork-info-in-alt-read-id "$bam_subset" |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py --window 300 --thres 0.5 --infer |\
    sed '1i\detectIndex mean_brdU start end' |\
    python filter_forks_autocorr_len_lessThanOrEqual.py --autoCorrThres 0.5 --lenThres 2 |\
    python get_model_fit.py --win 1000 \
        --fitFile "$fit_dir"/investigate_fit_params_"${tag}" --xMax 30 \
        --winTol 0.2
} > "$fit_dir"/investigate_fit_"${tag}"
