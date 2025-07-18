#!/bin/bash
#SBATCH --mem=120G
#SBATCH -c 1
#SBATCH -p ei-gpu
#SBATCH --gres=gpu:1
#SBATCH -J gup
#SBATCH --mail-type=END,FAIL
#SBATCH --time=72:00:00

set -e

# load git labels
source load_git_repo_labels.sh

# load guppy
source load_package.sh -"${guppyVersion:-"guppy_5_07"}"

guppy_basecaller --input_path $fast5Dir --save_path $fastQDir\
    --config $guppyConfigFile --recursive -q 4000 --disable_pings\
    --device cuda:0 --disable_qscore_filtering

# clean files up
cd $fastQDir
mkdir $guppyLogFileDir
mv ./*.log $guppyLogFileDir
cat ./*.fastq > tmp1234
rm -rf ./*.fastq
mv tmp1234 $fastQFile
