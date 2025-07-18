#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J getMnBPerFork
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# preamble
# --------

# get one (or one value per window over many windows) mean brdU value per fork called by forkSense
# usage sbatch <program_name.sh> $modBam $leftFork $rightFork $window_size_optional
# modBAM: .mod.bam file containing analogue modification probabilities
# leftFork: left fork file output by forksense
# rightFork: right fork file output by forksense
# window_size_optional: (optional) window along fork direction in number of thymidines if requested.
#                       defaults to one window per fork.

# Four columns are output to stdout with no header and separated by spaces
# index, mean_val, start, end.
# names are self-explanatory.

# stop execution if any command fails
set -e

# need at least 3 arguments
if [ $# -lt 3 ]
then
    >&2 echo "usage: sbatch <program_name.sh> \$modBam \$leftFork \$rightFork \$window_size_optional"
    >&2 echo "modBAM: .mod.bam file containing analogue modification probabilities"
    >&2 echo "leftFork: left fork file output by forksense"
    >&2 echo "rightFork: right fork file output by forksense"
    >&2 echo "window_size_optional: (optional) window along fork direction in number of thymidines if requested."
    >&2 echo "                      defaults to one window per fork."
    exit 1
fi

window_size=${4:-0}

# check that window size is greater than or equal to zero
if ! [ "$window_size" -ge 0 ]
then
    >&2 echo "window size must be greater than or equal to zero"
    exit 1
fi

# load git repo labels
source load_git_repo_labels.sh

# load python, miller
source load_package.sh -python -miller

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)

# set fork-location files and bam file
rtForks=$3
ltForks=$2
bam=$1

# if a fork file does not exist, then set it to /dev/null
if ! [ -f "$rtForks" ]
then
    rtForks=/dev/null
fi

if ! [ -f "$ltForks" ]
then
    ltForks=/dev/null
fi

# check that at least one of the fork files exists
if ! [ -f "$rtForks" ] && ! [ -f "$ltForks" ]
then
    >&2 echo "at least one input feature file must exist"
    exit 1
fi

# check that the bam file exists
if ! [ -f "$bam" ]
then
    >&2 echo "bam file does not exist"
    exit 1
fi

# check that the bam file is indexed
if ! [ -f "$bam".bai ]
then
    >&2 echo "bam file is not indexed"
    exit 1
fi

# make temporary directory
mkdir -p "$tmpDir"

# gather forks and split into 2000 forks per file
{
    < "$ltForks" awk '{print $1 " " $2 " " $3 " " $4 " " $4 "_" $5 "_" $6 "_" $7 "_" $8 "_L_" $2 "_" $3}';
    < "$rtForks" awk '{print $1 " " $2 " " $3 " " $4 " " $4 "_" $5 "_" $6 "_" $7 "_" $8 "_R_" $2 "_" $3}';
} | split -l 2000 - "$tmpDir"/forks

# do processing per split file
for splitFile in "$tmpDir"/forks??
do
        <"$splitFile" sed '1i\contig start end read_id alt_read_id' |\
        python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$bam" |\
        sed '1i\detectIndex\tposOnRef\tval' |\
        python get_mean_brdU_window.py --thres 0.5 --window "$window_size" --infer\
            > "$splitFile"_processed &
done

wait;

# set the script information
echo "# from commit ${COMMITSTR} generated at ${TIMENOW}";
echo "# script: $0";
echo "# arguments: $*";
mlr --csv --ifs ' ' --ofs ' ' --implicit-csv-header --headerless-csv-output \
  --skip-comments cat "$tmpDir"/forks??_processed;

rm -rf "$tmpDir"