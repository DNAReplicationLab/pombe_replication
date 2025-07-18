#!/bin/bash

#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J collatePauseToMsrPause
#SBATCH --time=5:59:59

# run the script using the command
# sbatch --dependency=afterok:$jid scriptName.sh $inputDir $jid $outputFile

# fail on error
set -e

# the script takes all output files with jid in them in the input directory,
# collates them into a file in the output directory, and deletes job files.

inputDir=${1:-}
jid=${2:-}
outputFile=${3:-}

# check that the input directory is specified
if [ -z "$inputDir" ]; then
    >&2 echo "Input directory not specified. Exiting."
    exit 1
fi

# check that the input directory exists
if [ ! -d "$inputDir" ]; then
    >&2 echo "Input directory does not exist. Exiting."
    exit 1
fi

# check that the jid is specified
if [ -z "$jid" ]; then
    >&2 echo "Job ID not specified. Exiting."
    exit 1
fi

# create the output file if it doesn't exist
if [ -z "$outputFile" ]; then
    >&2 echo "Output file not specified. Exiting."
    exit 1
else
    outputDir=$(dirname "$outputFile")
    mkdir -p "$outputDir"
fi

# find all files, if first file, then just send contents
# to output file, otherwise delete all leading comments
# and header and then send to file, then delete job file.

declare -i Q=0

# shellcheck disable=SC2154
for entry in "$inputDir"/*"$jid"*.out
do

    if [ $Q -eq 0 ]
    then
        cp "$entry" "$outputFile";
    else
        grep -v '^#' "$entry"  >> "$outputFile";
    fi

    rm "$entry";
    Q=1;

done
