#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 3
#SBATCH -p ei-medium
#SBATCH -J runPlotGetFirstLastWinMnBrduPerFork
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# calculate and plot histogram of mean BrdU density in first and last 1kb of forks

# usage
#------
# bash run_and_plot_get_first_last_win_mean_brdu_per_fork.sh <mod_bam> <forksense_dir> <op_file>
# <mod_bam> - path to modified bam file
# <forksense_dir> - path to forksense directory which contains origins, leftForks, rightForks etc.
# <op_file> - path to output file, will be created/overwritten

# outputs
# -------
# <op_file>, <op_file>.forkStart.png, <op_file>.forkEnd.png
# <op_file> is space-separated. For contents, see get_first_last_win_mean_brdu_per_fork.sh

# stop execution if any command fails
set -euxo pipefail

# load packages
source load_package.sh -miller -R

# load input variables
mod_bam=$1
forksense_dir=$2
op_file=$3
win_size=999        # in bp, the window size to calculate mean BrdU density in
n_forks_sample=5000 # number of forks to sample, keep it in 1000s, if you want more, ensure you understand
                    # this script, the script being called, and are prepared to increase memory and/or time,
                    # and inspect the output carefully.

# check if bam file exists
if ! [ -f "$mod_bam" ]
then
    >&2 echo "bam file does not exist"
    exit 1
fi

# check that the bam file is indexed
if ! [ -f "$mod_bam".bai ]
then
    >&2 echo "bam file is not indexed"
    exit 1
fi

# check if forksense directory exists
if ! [ -d "$forksense_dir" ]
then
    >&2 echo "forksense directory does not exist"
    exit 1
fi

# calculate mean BrdU density in first and last (approx) 1kb of forks
bash get_first_last_win_mean_brdu_per_fork.sh "$mod_bam" "$forksense_dir" $win_size $n_forks_sample > "$op_file"

# plot histogram of mean BrdU density in first window of forks
# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --csv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments rename 2,mean_brdu,5,type then\
  filter '$type == "first"' then\
  histogram -f mean_brdu --lo 0 --hi 1 --nbins 40 "$op_file" |\
 Rscript plotting_and_short_analyses/plot_histogram.R \
   mean_brdu_count "$op_file".forkStart.png "BrdU density in fork first kb" "Count" 0,1 0,auto

# plot histogram of mean BrdU density in last window of forks
# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --csv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments rename 2,mean_brdu,5,type then\
  filter '$type == "last"' then\
  histogram -f mean_brdu --lo 0 --hi 1 --nbins 40 "$op_file" |\
 Rscript plotting_and_short_analyses/plot_histogram.R \
   mean_brdu_count "$op_file".forkEnd.png "BrdU density in fork last kb" "Count" 0,1 0,auto
