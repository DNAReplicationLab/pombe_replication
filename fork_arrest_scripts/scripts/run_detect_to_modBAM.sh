#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J detectToModBAM
#SBATCH --mail-type=END,FAIL
#SBATCH --time 47:59:59

# preamble
# --------

# stop execution if any command fails
set -e

helpFunction()
{
   echo "Convert detect files to modBAM and sort and index them"
   echo "Usage: $0 -f detectFile"
   echo -e "\t-f /path/to/detectFile.detect"
   exit 1 # Exit script after printing help
}

while getopts "f:" opt
do
   case "$opt" in
      f ) detectFile="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case param doesn't exist
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$detectFile" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction;
fi

# load git repo labels
source load_git_repo_labels.sh

# load python, samtools
source load_package.sh -python -samtools

modBAMPrefix=$detectFile.mod

# set the script information
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# set info about BrdU modification
infoMod="# using generic sam format thymidine modification tag T to mean BrdU"

# perform conversion
< "$detectFile" \
    sed "1i${infoStr}\n${infoMod}" |\
    python DNAscentTools/convert_detect_to_modBAM.py \
      --op "$modBAMPrefix".bam --tag T

# index and sort files
samtools sort -o "$modBAMPrefix".sorted.bam "$modBAMPrefix".bam
samtools index "$modBAMPrefix".sorted.bam
