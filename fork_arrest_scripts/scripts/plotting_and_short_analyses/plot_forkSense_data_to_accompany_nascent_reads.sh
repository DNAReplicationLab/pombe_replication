#!/bin/bash

#SBATCH --mem=5G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotFrkForNasReads
#SBATCH --mail-type=END,FAIL
#SBATCH --time 4:59:59

# goal of program
# =================
# Obtain read ids of plots in folder and produce one left and right forksense raw data plot per read.
# Prerequisite: Run main pipeline of fast5 -> dnascent before this.

# program sample usage and input parameter description
# ====================================================
# sbatch <program.sh> mod_bam_left mod_bam_right op_dir
# mod_bam_left, mod_bam_right : BAM file containing modification data. will have been produced by pipeline.
# op_dir: directory where the n plots go

# stop execution if any command fails
set -e

# load variables from the command line
mod_bam_left=$1;
mod_bam_right=$2;
op_dir=$3;

# obtain read ids in the directory and plot raw data of left and right fork probabilities
# not all read ids have forkSense data associated with them, so some of these jobs will fail and it's not a problem.
find "$op_dir"/plot*.png | awk -F"_" '{print $NF}' |\
    awk -F. '{print "sbatch --mail-user= -J frkNas limited_plot_read.sh '"$mod_bam_left"' " $1 " 0 '"$op_dir"' forkSense_left"}' |\
    bash

find "$op_dir"/plot*.png | awk -F"_" '{print $NF}'  |\
    awk -F. '{print "sbatch --mail-user= -J frkNas limited_plot_read.sh '"$mod_bam_right"' " $1 " 0 '"$op_dir"' forkSense_right"}' |\
    bash

# launch a dummy job so that we wait till plot jobs have finished
sbatch --dependency=singleton --job-name=frkNas --mail-user= --wait -p ei-short --wrap="sleep 10;"