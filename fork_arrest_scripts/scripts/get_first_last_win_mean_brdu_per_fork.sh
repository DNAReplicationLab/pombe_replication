#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 3
#SBATCH -p ei-medium
#SBATCH -J getFirstLastWinMnBrduPerFork
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# ----
# get the mean BrdU incorporation in the first and last windows of each fork

# usage
# -----
# sbatch <program_name.sh> $modBam $forkSense_dir $window_size_optional $sample_forks_optional
# modBAM: .mod.bam file containing analogue modification probabilities
# forkSense_dir: directory containing the forksense files with standard names and in standard format
# window_size_optional: (optional) measure over this window size (bp) at start and end of fork, default 1000.
# sample_forks_optional: (optional) sample this many forks at random (to save compute time), default = 1000.
#                        if you want a high number, then adjust number of cores and inspect output.

# output
# ------
# Four columns are output to stdout with no header and separated by spaces
# index, mean_val, start, end, label.
# names are self-explanatory, except for label, which is either "first" or "last" depending on whether the
# window is at the start or end of the fork.

# stop execution if any command fails
set -e

# need at least 2 arguments
if [ $# -lt 2 ]
then
    >&2 echo "usage: sbatch <program_name.sh> \$modBam \$forkSense_dir \$window_size_optional \$sample_forks_optional"
    >&2 echo "modBAM: .mod.bam file containing analogue modification probabilities"
    >&2 echo "forkSense_dir: directory containing the forksense files with standard names and in standard format"
    >&2 echo "window_size_optional: (optional) measure over this window size (bp) at fork start and end, default 1000"
    >&2 echo "sample_forks_optional: (optional) sample this many forks at random (to save compute time), default = 1000."
    >&2 echo "                        if you want a high number, then adjust number of cores and inspect output."
    exit 1
fi

sample_forks_optional=${4:-1000}
window_size=${3:-1000}
min_fork_len_kb=2
min_align_len_kb=30

# check that window size is an integer
if ! [[ "$window_size" =~ ^[0-9]+$ ]]
then
    >&2 echo "window size must be an integer."
    exit 1
fi

# check that window size is greater than zero
if ! [ "$window_size" -gt 0 ]
then
    >&2 echo "window size must be greater than zero"
    exit 1
fi

# check that window size is less than half the minimum fork length
if [ "$(echo "$window_size >= $min_fork_len_kb * 1000 / 2" | bc -l)" -eq 1 ]
then
    >&2 echo "window size must be less than twice the minimum fork length"
    exit 1
fi

# check that sample_forks_optional is an integer
if ! [[ "$sample_forks_optional" =~ ^[0-9]+$ ]]
then
    >&2 echo "sample_forks_optional must be an integer."
    exit 1
fi

# check that sample_forks_optional is greater than zero
if ! [ "$sample_forks_optional" -gt 0 ]
then
    >&2 echo "sample_forks_optional must be greater than zero"
    exit 1
fi

# load git repo labels
source load_git_repo_labels.sh

# load python, miller, jq
source load_package.sh -python -miller -jq

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)

# set fork-location and bam file
bam=$1
forkSense_dir=$2

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

# filter forkSense files
cd specialized_analyses;
bash filter_forksense_files_on_length_and_noChrM.sh "$forkSense_dir" "$min_fork_len_kb" \
  "$min_align_len_kb" "$tmpDir" > "$tmpDir"/filter_forksense_files_on_length_and_noChrM.json
cd ..

# get list of subset left and right forks
ltForks=$(< "$tmpDir"/filter_forksense_files_on_length_and_noChrM.json jq -r '.left_fork_subset_file')
rtForks=$(< "$tmpDir"/filter_forksense_files_on_length_and_noChrM.json jq -r '.right_fork_subset_file')

# gather forks, sample them, and split into 2000 entries per file
n_sample_windows=$((sample_forks_optional * 2))
{
    < "$ltForks" awk -v OFS=' ' -v w="$window_size" '{print $1, $3 - w, $3, $4, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_L_" $2 "_" $3 ":first:"}';
    < "$ltForks" awk -v OFS=' ' -v w="$window_size" '{print $1, $2, $2 + w, $4, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_L_" $2 "_" $3 ":last:"}';
    < "$rtForks" awk -v OFS=' ' -v w="$window_size" '{print $1, $3 - w, $3, $4, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_R_" $2 "_" $3 ":last:"}';
    < "$rtForks" awk -v OFS=' ' -v w="$window_size" '{print $1, $2, $2 + w, $4, $4 "_" $5 "_" $6 "_" $7 "_" $8 "_R_" $2 "_" $3 ":first:"}';
} | shuf | head -n "$n_sample_windows" | split -l 2000 - "$tmpDir"/forks

# do processing per split file
for splitFile in "$tmpDir"/forks??
do
       < "$splitFile" sed '1i\contig start end read_id alt_read_id' |\
        python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$bam" |\
        sed '1i\detectIndex\tposOnRef\tval' |\
        python get_mean_brdU_window.py --thres 0.5 > "$splitFile"_processed &
done

wait;

# set the script information and send output to stdout
echo "# from commit ${COMMITSTR} generated at ${TIMENOW}";
echo "# script: $0";
echo "# arguments: $*";
mlr --csv --ifs ' ' --ofs ' ' --implicit-csv-header --headerless-csv-output \
  --skip-comments cat "$tmpDir"/forks??_processed | awk -F: '{print $1  $3 " " $2}'

# remove temporary directory
rm -rf "$tmpDir"