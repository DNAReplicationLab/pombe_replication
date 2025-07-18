#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 16
#SBATCH -p ei-medium
#SBATCH -J nasindx
#SBATCH --mail-type=END,FAIL
#SBATCH --time=09:59:59

set -e

# load miller
source load_package.sh -miller

if [ "${useDnascentVersion4:-0}" -eq 1 ]; then
    # load dnascent v4
    source load_package.sh -dnascent_v4
else
    # load dnascent
    source load_package.sh -dnascent
fi

# make a temporary file, keeping only the column we are interested in
# shellcheck disable=SC1010
mlr --itsv --otsv rename filename,filename_fast5 \
  then cut -f filename_fast5,read_id "$guppySequencingSummaryFile" > "$guppySequencingSummaryFile"_rmCol

# perform indexing
DNAscent index -f "$fast5Dir" -o "$dnascentIndexFile" -s "$guppySequencingSummaryFile"_rmCol

# remove the temporary file
rm -rf "$guppySequencingSummaryFile"_rmCol