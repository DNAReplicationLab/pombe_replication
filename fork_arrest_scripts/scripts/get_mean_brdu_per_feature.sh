#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J getMnBPerFeature
#SBATCH --mail-type=END,FAIL
#SBATCH --time 23:59:59
#SBATCH --constraint=""

# preamble
# --------

# get one (or one value per window over many windows) mean brdU value per feature called by forkSense
# usage sbatch <program_name.sh> $modBam $featureFile $window_size_optional
# modBAM: .mod.bam file containing analogue modification probabilities
# featureFile: feature file output by forksense
# window_size_optional: (optional) window data in number of thymidines if requested.
#                       defaults to one window per feature.

# Four columns are output to stdout with no header and separated by spaces
# index, mean_val, start, end.
# names are self-explanatory.

# stop execution if any command fails
set -e

# need at least two arguments
if [ $# -lt 2 ]
then
    >&2 echo "usage: sbatch <program_name.sh> \$modBam \$featureFile \$window_size_optional"
    >&2 echo "modBAM: .mod.bam file containing analogue modification probabilities"
    >&2 echo "featureFile: feature file output by forksense"
    >&2 echo "window_size_optional: (optional) window data in number of thymidines if requested."
    >&2 echo "                      defaults to one window per feature."
    exit 1
fi

# use get_mean_brdu_per_fork.sh to get mean brdU per feature, but replace _R_ with _F_ in the output
# because the feature file used as an input here may not necessarily be a right fork file
bash get_mean_brdu_per_fork.sh "$1" /dev/null "$2" "${3:-0}" | sed 's/_R_/_F_/g'