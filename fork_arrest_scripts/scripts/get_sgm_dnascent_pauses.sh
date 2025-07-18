#!/bin/bash
#SBATCH --mem=6G
#SBATCH -c 4
#SBATCH -p ei-medium
#SBATCH -J getSgmDnascentPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --time 4:59:59

# goal
# -----
# Start from dnascent data and perform pause detection and processing

# usage
# -----
# sbatch get_sgm_dnascent_pauses.sh json_file

# json_file: json file with parameters for the pipeline
# the json_file should look as follows (the fields can be in any order and any number of additional fields are ok):
# {
#   "left_forks": "/path/to/left_forks_file.forksense",
#   "right_forks": "/path/to/right_forks_file.forksense",
#   "mod_bam": "/path/to/mod_bam_file.bam",
#   "output_dir": "/path/to/output/directory",
#   "output_file_prefix": "output_file_prefix",
#   "low": 0.1,
#   "high": 0.7,
#   "width": 3400,
#   "mod_bam_left": "/path/to/mod_bam_left.bam",
#   "mod_bam_right": "/path/to/mod_bam_right.bam",
#   "pipeline_go_get_model_pauses": 1,
#   "pipeline_go_perform_further_analysis": 1
# }

# explanation of json fields:
# ltForks: left forks file, in forkSense format
# rtForks: right forks file, in forkSense format
# modBam: modified bam file, obtained by conversion from dnascent detect output
# outputDir: output directory
# outputFilePrefix: prefix for output files i.e. output files are <outputDir>/<outputFilePrefix>.<some_extension>
# low: lowest level of sigmoid
# high: highest level of sigmoid
# width: width of sigmoid
# modBamLeft: modified bam file containing fork probabilities per position per read for left forks
# modBamRight: modified bam file containing fork probabilities per position per read for right forks
# pipelineGoGetModelPauses: (default 0) set to 1 if you want to run the pause getting script
# pipelineGoUseTheoreticalFilter: deprecated
# pipelineGoPerformFurtherAnalysis: (default 0) set to 1 if you want to run the further analysis

# fail on error
set -e

# load git labels
source load_git_repo_labels.sh

# load jquery
source load_package.sh -jq

# load some config information
source config.sh

# exit if email information doesnt exist
# shellcheck disable=SC2154
if [ "${config[email]}" == "" ]; then
    echo "email not set in config.sh"
    exit 1
fi

# check that input arguments have been supplied
if [ "$#" -ne 1 ]; then
    >&2 echo "Illegal number of parameters"
    >&2 echo "usage: bash get_sgm_dnascent_pauses.sh <json_file>"
    >&2 echo "json_file: json file with parameters for the pipeline"
    >&2 echo "see the header of this script for an example of the json file"
    exit 1
fi

# check that the json file exists
if [ ! -f "$1" ]; then
    echo "json file does not exist" >&2; exit 1
fi

# set input data
ltForks=$(jq -r '.left_forks' "$1")
rtForks=$(jq -r '.right_forks' "$1")
bam=$(jq -r '.mod_bam' "$1")
outputDir=$(jq -r '.output_dir' "$1")
outputFilePrefix=$(jq -r '.output_file_prefix' "$1")
pauseFile="$outputDir"/"$outputFilePrefix".pauses
tempPauseFile="$pauseFile"_"$(openssl rand -hex 6)"
low=$(jq -r '.low' "$1")
high=$(jq -r '.high' "$1")
width=$(jq -r '.width' "$1")
bam_left=$(jq -r '.mod_bam_left' "$1")
bam_right=$(jq -r '.mod_bam_right' "$1")

# check that the input files exist
if [ ! -f "$ltForks" ] || [ ! -f "$rtForks" ]; then
    echo "fork file(s) do not exist" >&2; exit 1
fi

if [ ! -f "$bam" ] || [ ! -f "$bam_left" ] || [ ! -f "$bam_right" ]; then
    echo "bam file(s) do not exist" >&2; exit 1
fi

# check that the bam index files exist
if [ ! -f "$bam".bai ] || [ ! -f "$bam_left".bai ] || [ ! -f "$bam_right".bai ]; then
    echo "bam index file(s) do not exist" >&2; exit 1
fi

# check that low, high and width are numbers
re='^[0-9]+([.][0-9]+)?$'

if ! [[ $low =~ $re ]] || ! [[ $high =~ $re ]] || ! [[ $width =~ $re ]]; then
   echo "error: either low, high, or width, or all are not numbers" >&2; exit 1
fi

# check that the output directory exists, otherwise create it
if [ ! -d "$outputDir" ]; then
    mkdir -p "$outputDir"
fi

# first, subset the bam files to retain only the reads that are in the fork files
echo "Subsetting bam files to retain only the reads that are in the fork files"
read_ids="$pauseFile".read_ids
tmp_file=$(mktemp -p "$outputDir")
cat "$ltForks" "$rtForks" | awk '{print $4}' | sort | uniq > "$read_ids"
bash subset_four_bam_files.sh /dev/null "$bam" "$bam_left" "$bam_right" "$read_ids" "$outputDir" > "$tmp_file"

bam=$(jq -r '.mod_bam_subset' "$tmp_file")
bam_left=$(jq -r '.mod_bam_left_subset' "$tmp_file")
bam_right=$(jq -r '.mod_bam_right_subset' "$tmp_file")
rm "$tmp_file"

# do a dummy job to get an initial job id
launchJob=$(sbatch --wrap="sleep 1" -p ei-short --time=00:30:10\
    -o /dev/null -e /dev/null)
jid=${launchJob##* }

# set the variable below to one before the section of the pipeline
# that needs to be run. set it to zero before whichever part of the
# pipeline that needn't be run.
PIPELINE_GO=$(jq -r '.pipeline_go_get_model_pauses' "$1")
if [ "$PIPELINE_GO" -ne 0 ] && [ "$PIPELINE_GO" -ne 1 ]; then
    echo "pipeline_go_get_model_pauses should be 0 or 1, so setting it to 0" >&2
    PIPELINE_GO=0;
fi

# find pauses using our cut and align procedure
if [ "$PIPELINE_GO" -eq 1 ];
then

    # count total number of forks by adding the number of lines in the left and right forks files
    n_forks=$(wc -l < "$ltForks")
    n_forks=$((n_forks + $(wc -l < "$rtForks")))

    # if there are no forks, raise an error
    if [ "$n_forks" -eq 0 ]; then
        echo "no forks found in the input files" >&2;
        exit 1
    fi

    # round up to the nearest multiple of 1000
    n_forks_in_multiples_of_1000=$(( (n_forks + 999) / 1000 ))
    n_forks_in_multiples_of_1000=$((n_forks_in_multiples_of_1000 - 1)) # job array numbering starts from 0

    launchJob=$(sbatch --dependency=afterok:"$jid"\
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       --array=0-${n_forks_in_multiples_of_1000} \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_get_model_pauses_raw_sgm_process.sh "$ltForks" "$rtForks" "$bam" "$tempPauseFile"\
                            "$low" "$high" "$width" ;)
    jid=${launchJob##* }

    # the previous step produces a number of output files.
    # we have to stitch them together.
    # Important note: this stage runs even if the previous step has failed.
    #                 This is because the previous stage runs many jobs in parallel, and we cannot afford to stop
    #                 the whole pipeline if one of the jobs fails.

    launchJob=$(sbatch --dependency=afterany:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       stitch_many_pause_files_together.sh "$tempPauseFile" "$pauseFile" ;)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(jq -r '.pipeline_go_use_theoretical_filter // "NA"' "$1")
if [ "$PIPELINE_GO" != "NA" ]; then
    echo "pipeline_go_use_theoretical_filter is deprecated" >&2
fi

PIPELINE_GO=$(jq -r '.pipeline_go_perform_further_analysis' "$1")
if [ "$PIPELINE_GO" -ne 0 ] && [ "$PIPELINE_GO" -ne 1 ]; then
    echo "pipeline_go_perform_further_analysis should be 0 or 1, so setting it to 0" >&2
    PIPELINE_GO=0;
fi

# perform analysis of our pause results
if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_parallel_run_process_model_pause_results_analysis.sh "$pauseFile"\
                        "$pauseFile"_after_analysis "$bam" "$bam_left" "$bam_right" ;)
    jid=${launchJob##* }
fi
