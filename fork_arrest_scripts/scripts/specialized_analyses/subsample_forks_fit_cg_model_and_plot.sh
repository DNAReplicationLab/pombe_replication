#!/bin/bash

#SBATCH --mem=4G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J fitCGModel
#SBATCH --mail-type=END,FAIL
#SBATCH --time 3:59:59

# goal
# --------
# filter forksense left and right fork files, and fit model
# to windowed data from N forks selected at random

# description
# ------------
# after our main pipeline is run, we get a mod bam file and some fork sense files.
# we are going to filter the left and right fork files according to some criteria,
# say "preserve only forks with lengths above 10 kb found on reads with alignments above 30 kb and avoid chrM."
# NOTE: the criteria is not mentioned here as it could change.
# then, we will select N/2 forks at random from the left and right forksense files,
# (where N is a user supplied number), and fit our coarse grained model to it.
# (by coarse grained, we mean we are going to window the fork data on genomic-coordinate windows
#  and fit a model to it, as opposed to fitting a model to the un-windowed data).

# usage and input
# ---------------
# bash <program_name.sh> $mod_bam $forkSense_directory $op_dir $tag $n_forks
# NOTE: use sbatch for running this script if you are on a cluster and want to sample many forks.
# * inputs are (1) a mod bam file with analogue modification probabilities,
# *            (2) the fork sense directory with the files with left forks, right forks etc.
# *            (3) output directory where fit results and plots will go
# *            (4) a short string used as a suffix to tag results
# *            (5) number of forks used in model fitting. use a number in the 100s - max 1000 or -1 (read below).
# *                for numbers above 1000, try a few things like running the script as is, changing memory etc.
# *                if these don't work, then a new approach / modifications may be needed.
# *                If you use -1, then the program expects a suitably named leftForks and rightForks file to be
# *                present in the output directory, and will use those files for fitting. This is useful if you want
# *                to re-run the program again in the future without resampling, or if you have a specific set of forks
# *                you want to fit the model to.
# *                If you use -1, then the program will not randomly subsample forks, so the fork sense
# *                directory parameter is unused.
# * the left and right forks file must have the standard names leftForks_DNAscent_forkSense.bed,
#   rightForks_DNAscent_forkSense.bed, which are the names given by fork sense.
#   They must have the standard fork sense format i.e. space separated, no headers,
#   with columns contig,start_fork,end_fork,read_id,contig_again,start_align,end_align,orn.

# read command line inputs
if [ "$#" -ne 5 ]; then
    echo "Usage: bash <program_name.sh> mod_bam forkSense_directory op_dir tag n_forks"
    exit;
fi

# set directories
pwd=$(pwd)
plotting_dir=../plotting_and_short_analyses
config_dir=..
script_dir=..

# load configuration variables and packages
cd "$config_dir" || exit;
source config.sh
source load_package.sh -miller -R
cd "$pwd" || exit;

# accept inputs
mod_bam=$1

mkdir -p "$3";
fit_dir=$(cd "$3" || exit; pwd)

tag=$4
if [ "$5" -eq -1 ]; then
  echo "Using the left and right forks files in the output directory for fitting"
  n_forks=-1
  forksense_dir=""
else
  n_forks=$(echo "$5"/2 | bc)
  forksense_dir=$(cd "$2" || exit; pwd)
fi

# perform filtration and sampling on forksense files if n_forks is not -1
for dirn in "left" "right"
do
  output_file="$fit_dir"/"$dirn"Forks_DNAscent_forkSense.forkLen.10kb.alignLen.30kb.noChrM.subset."$tag".bed

  if [ "$n_forks" -eq -1 ]; then
    if [ ! -f "$output_file" ]; then
      echo "The file $output_file does not exist. So you cannot use -1 for n_forks. Exiting."
      exit
    else
      continue
    fi
  fi

  input_file="$forksense_dir"/"$dirn"Forks_DNAscent_forkSense.bed

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --headerless-csv-output \
    label contig,start_fork,end_fork,read_id,contig_again,start_align,end_align,orn \
    then put '$fork_len = $end_fork - $start_fork; $align_len = $end_align - $start_align' \
    then filter '$fork_len >= 10000 && $align_len >= 30000 && $contig != "chrM"' \
    then cut -x -f fork_len,align_len \
    then shuffle "$input_file" > "$output_file".tmp

  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --headerless-csv-output head -n "$n_forks" \
    "$output_file".tmp > "$output_file"

  echo "After filtering $dirn forks, we found $(wc -l < "$output_file".tmp) forks and sampled $n_forks"

  rm "$output_file".tmp

done

# window data and fit model to forks
cd "$script_dir" || exit;
bash run_get_model_fit.sh "$mod_bam"\
  "$fit_dir"/leftForks_DNAscent_forkSense.forkLen.10kb.alignLen.30kb.noChrM.subset."$tag".bed\
  "$fit_dir"/rightForks_DNAscent_forkSense.forkLen.10kb.alignLen.30kb.noChrM.subset."$tag".bed\
  "$fit_dir" "$tag";
cd "$pwd" || exit;

# plot the fits
cd "$plotting_dir" || exit;
Rscript plot_model_fit_fork_data.R "${fit_dir}"/investigate_fit_"${tag}" "${fit_dir}"/investigate_fit_"${tag}".png