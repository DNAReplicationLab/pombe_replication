#!/bin/bash
#SBATCH --mem=120G
#SBATCH -c 16
#SBATCH -p ei-long
#SBATCH -J dnascDet
#SBATCH --mail-type=END,FAIL
#SBATCH --time=199:00:00

set -e

if [ "${useDnascentVersion4:-0}" -eq 1 ]; then
    # load dnascent v4
    source load_package.sh -dnascent_v4
else
    # load dnascent
    source load_package.sh -dnascent
fi

numCPUThreads=30

# need HDF5 plugin for newer fast5 formats
# shellcheck disable=SC2154
export HDF5_PLUGIN_PATH=$HDF5PluginPath

# shellcheck disable=SC2154
DNAscent detect -b "$minimap2SortedBamFile" -r "$refGenomeFasta" \
 -i "$dnascentIndexFile" -o "$dnascentDetectFile" -t $numCPUThreads \
 -q "$minQual" -l "$minLen"
