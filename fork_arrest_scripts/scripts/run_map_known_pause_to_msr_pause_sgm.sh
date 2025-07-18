#!/bin/bash

#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J mapPauseToMsrPause
#SBATCH -o /hpc-home/thiyagar/fork_arrest/scratch/tmp/mapPauseToMsrPause.slurm.%N.%A.%a.out
#SBATCH -e /hpc-home/thiyagar/fork_arrest/scratch/tmp/mapPauseToMsrPause.slurm.%N.%A.%a.err
#SBATCH --constraint=""
#SBATCH --time=2:59:59
#SBATCH --array=1-1000

# script is an array job. after completion, a collation
# job has to be run to gather all results into one file.
# collation script is separate.

# load python
source load_package.sh -python

# print requisite fork msr to known pause maps to output
python map_known_pause_to_msr_pause_sgm.py -n 100 -high 0.65 -low 0.1 -numeric -iter 20
