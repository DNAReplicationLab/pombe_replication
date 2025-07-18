#!/bin/bash

#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J forkSenseToModBAM
#SBATCH --mail-type=END,FAIL
#SBATCH --time 47:59:59
#SBATCH --constraint=""

# preamble
# --------

# stop execution if any command fails
set -e

helpFunction()
{
   echo "Convert forkSense files to modBAM and sort and index them"
   echo "Usage: $0 -f forkSenseFile -r referenceFasta -d left"
   echo -e "\t-f /path/to/forkSenseFile.forkSense"
   echo -e "\t-r /path/to/reference.fasta"
   echo -e "\t-d left or -d right. convert left or right fork probabilities"
   echo "Output file is /path/to/forkSenseFile.forkSense.mod.left.sorted.bam or same name with right instead of left"
   exit 1 # Exit script after printing help
}

while getopts "f:r:d:" opt
do
   case "$opt" in
      f ) forkSenseFile=${OPTARG:?word} ;;
      r ) reference=${OPTARG:?word} ;;
      d ) direction=${OPTARG:?word} ;;
      ? ) helpFunction ;; # Print helpFunction in case param doesn't exist
   esac
done

if ! { [ "$direction" == "left" ] || [ "$direction" == "right" ]; };
then
  echo 'Direction should be left or right';
  helpFunction ;
  exit 1;
fi

# check that the forkSense file and the reference file exist
if [ ! -f "$forkSenseFile" ];
then
  echo "forkSense file $forkSenseFile does not exist";
  helpFunction ;
  exit 1;
fi

if [ ! -f "$reference" ];
then
  echo "reference file $reference does not exist";
  helpFunction ;
  exit 1;
fi

# load git repo labels
source load_git_repo_labels.sh

# load python, samtools
source load_package.sh -python -samtools

# set output filename
modBAMPrefix=${forkSenseFile}.mod.${direction}

# set the script information
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

# set info about BrdU modification
infoMod="# using generic sam format thymidine modification tag T"

# info about forkSense
infoFs="# as this is a ($direction) forkSense file, modification probabilities must be taken to mean fork probabilities"

if [ "$direction" == "left" ];
then
  convert_option=;
else
  convert_option="--switchCols2and3";
fi

# perform conversion
< "$forkSenseFile" \
    sed "1i${infoStr}\n${infoMod}\n${infoFs}" |\
    python DNAscentTools/convert_detect_to_modBAM.py \
      --op "$modBAMPrefix".bam --tag T --fasta "$reference" \
      $convert_option
# use the switch option if necessary

# index and sort files
samtools sort -o "$modBAMPrefix".sorted.bam "$modBAMPrefix".bam
samtools index "$modBAMPrefix".sorted.bam
