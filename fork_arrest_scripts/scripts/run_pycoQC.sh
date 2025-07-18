#!/bin/bash

#SBATCH --mem=20G
#SBATCH -p ei-medium
#SBATCH -J pycoQC
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

if [ "$#" -ne 4 ]; then
    echo -n "Usage: bash scriptName.sh inputSortedBam inputGuppySummaryFile "
    echo "outputPycoQCDir tempDirPath"
    echo -e "\t* Can also use sbatch in place of bash."
    echo -e "\t* outputPycoQCDir, tmpDir will be created if they don't exist."
    exit 2;
fi

minimap2SortedBamFile=$1
summaryFile=$2
pycoQCAnalysisDir=$3
tmpDir=$4

# check that the sequencing summary file and the bam file exist
if [ ! -f "$summaryFile" ]; then
    echo "sequencing summary file $summaryFile does not exist" >&2
    exit 1
fi
if [ ! -f "$minimap2SortedBamFile" ]; then
    echo "pycoQC bam file $minimap2SortedBamFile does not exist" >&2
    exit 1
fi

# load pycoQC
source load_package.sh -pycoQC

# make temp directory
mkdir -p "$tmpDir"

# make analysis directory
mkdir -p "$pycoQCAnalysisDir"

# run pycoQC
pycoQC \
    -f "$summaryFile" \
    -a "$minimap2SortedBamFile" \
    -o "$pycoQCAnalysisDir"/analysis.html \
    -j "$pycoQCAnalysisDir"/analysis.json \
    --quiet