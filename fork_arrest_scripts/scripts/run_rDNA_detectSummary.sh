#!/bin/bash

#SBATCH --mem-per-cpu=200G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J run_rDNA_detectSummary
#SBATCH --mail-type=END,FAIL
#SBATCH --time=47:59:59

# goal
# ====
# Find pauses in the rDNA region of the genome and perform some analysis and plotting on them.

# usage
# =====
# sbatch run_rDNA_detectSummary.sh $mod_bam $output_directory
# mod_bam: BAM file with modified base information. A BAM index file must be present in the same directory.
# output_directory: directory to store output files

# output
# ======
# a bunch of files are sent to the output directory

# stop execution if any step fails
set -e

# check if correct number of arguments were provided
if [ "$#" -lt 2 ]; then
  >&2 echo "Incorrect number of parameters"
  >&2 echo "Usage: sbatch run_rDNA_detectSummary.sh <mod_bam> <output_directory> <several_optional_parameters>"
  >&2 echo "mod_bam: BAM file with modified base information. A BAM index file must be present in the same directory."
  >&2 echo "output_directory: directory to store output files"
  >&2 echo "Following parameters are optional."
  >&2 echo "NOTE: Remember that if you want to specify an optional parameter, you must specify all the parameters before it in the order they are listed here."
  >&2 echo "      You can specify \"\" as a parameter(s) to bypass this requirement."
  >&2 echo "      e.g. sbatch run_rDNA_detectSummary.sh \$mod_bam \$output_directory \"\" \$contig_name "
  >&2 echo "      this will use the default value for ref_fasta"
  >&2 echo "ref_fasta: reference fasta file (default: config[dataDir]/references/sc_rdna_22_repeats_rc.fasta)"
  >&2 echo "contig_name: name of the contig in the reference fasta file (default: rDNA_22_repeat)"
  >&2 echo "num_repeat_units: number of repeat units in the reference fasta file (default: 22)"
  >&2 echo "bed_file_with_features: bed file with features (default: config[dataDir]/references/rDNA_reference_regions_on_sc_rdna_1_repeat.bed)"
  >&2 echo "base: base to use for pause aggregation into bedgraphs (default: N, which means all bases)"
  >&2 echo "number: number of bases to use for pause aggregation (default: 13)"
  >&2 echo "pos_boundary: position boundary to use for pause aggregation (default: 8855)"
  >&2 echo "pauseLocThresKb: pause-location-from-end threshold in kb to mark discards (default: 5)"
  >&2 echo "n: the s.d. prefactor to filter pauses by, must be an integer (default: 3)"
  exit 1;
fi

# set input and output directories/files
modBAM=$1
outputDir=$2

# check if input files exist
if [ ! -f "$modBAM" ]; then
  >&2 echo "Input file $modBAM does not exist"
  exit 1;
fi

# check that the bam index file exists
if [ ! -f "${modBAM}.bai" ]; then
  >&2 echo "BAM index file ${modBAM}.bai does not exist"
  exit 1;
fi

# check if output directory exists, if not create it
if [ ! -d "$outputDir" ]; then
  mkdir -p "$outputDir"
fi

# get config directory info
source config.sh

# check that config dataDir exists
if [ -z "${config[dataDir]}" ] || [ ! -d "${config[dataDir]}" ]; then
  >&2 echo "config[dataDir] is unset or doesn't exist";
  exit 1;
fi

# set other parameters
ref_fasta=${3:-"${config[dataDir]}"/references/sc_rdna_22_repeats_rc.fasta}
contig_name=${4:-rDNA_22_repeat}
num_repeat_units=${5:-22}
bed_file_with_features=${6:-"${config[dataDir]}"/references/rDNA_reference_regions_on_sc_rdna_1_repeat.bed}
base=${7:-N}
number=${8:-13}
pos_boundary=${9:-8855}
pauseLocThresKb=${10:-5}
n=${11:-3}

# check that ref_fasta and bed_file_with_features exist
if [ ! -f "$ref_fasta" ]; then
    >&2 echo "Reference fasta file $ref_fasta does not exist"
    exit 1;
fi

if [ ! -f "$bed_file_with_features" ]; then
    >&2 echo "Bed file with features $bed_file_with_features does not exist"
    exit 1;
fi

# check that n is an integer
if ! [[ "$n" =~ ^[0-9]+$ ]]; then
  >&2 echo "n must be an integer"
  exit 1;
fi

# load python, samtools
source load_package.sh -python -samtools -miller

# load git repo labels
source load_git_repo_labels.sh

{

  # print some info about the run
  echo "# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>";
  echo "# using $modBAM";
  echo "# slurm job name: ${SLURM_JOB_NAME:-NA}";

  # find pauses
  samtools view -h -e "rlen>=30000" "$modBAM" | python convert_modBAM_to_detect.py --fasta "$ref_fasta" |\
    python rDNA_detectSummary.py "$outputDir" "$n";

} > "$outputDir"/job_summary.txt

# process the nascent files
python rDNA_nascentStepSummary.py "$outputDir"/nascent.txt > "$outputDir"/nascentStepList.bed;
python rDNA_nascentStepSummary.py "$outputDir"/"$n"SD_step.txt | awk -v n="$n" '{if($5 > n){print $0}}' \
  > "$outputDir"/"$n"SD_step.split_by_step_and_filter_"$n"SD.bed;

# take the step data and convert it into a pause file for subsequent analysis
bash run_intersect_exact_bed_files.sh "$outputDir"/"$n"SD_step.split_by_step_and_filter_"$n"SD.bed \
    "$outputDir"/nascentAllStepListRaw.bed "$outputDir"/"$n"SD_step.split_by_step_and_filter_"$n"SD.addStepData.bed;

{

  # print some info about the run
  echo "# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>";
  echo "# using $outputDir/${n}SD_step.split_by_step_and_filter_${n}SD.bed";
  echo "# using $outputDir/nascentAllStepListRaw.bed";
  echo "# slurm job name: ${SLURM_JOB_NAME:-NA}";

  < "$outputDir"/"$n"SD_step.split_by_step_and_filter_"$n"SD.addStepData.bed python rDNA_process_calc_steps_into_pause_file.py

} > "$outputDir"/pauseList;

# make fake forksense calls for use with our plotting program
mkdir -p "$outputDir"/fakeForksenseFormat
< "$outputDir"/pauseList python convert_pause_file_to_forkSense_calls.py --outputDir "$outputDir"/fakeForksenseFormat

# make zero-file-size origin and termination forkSense calls
touch "$outputDir"/fakeForksenseFormat/origins_DNAscent_forkSense.bed
touch "$outputDir"/fakeForksenseFormat/terminations_DNAscent_forkSense.bed

# prepare temporary files
temp_file_1=$(mktemp)
temp_file_2=$(mktemp)

# prepare the pause file for plotting
bash rDNA_prepare_pauseFile_for_plotting.sh "$outputDir"/pauseList "$temp_file_1"

# get minimum and maximum analogue density per window
bash run_get_min_max_analogue_sliding_window_pauseFile.sh "$temp_file_1" "$modBAM" > "$temp_file_2"

# 1. mark pauses for discard where pause locations are too close to ends
# 2. mark pauses for discard where at least one sliding window has a density lower or higher than a threshold
# shellcheck disable=SC2016
python mark_forks_pauses_at_align_ends.py --pauseLocThresKb "$pauseLocThresKb" "$temp_file_2" |\
  mlr --tsv put '$keep_min_300T_sliding_window_along_fork_gt_0pt01 = capitalize(string($min_300T_sliding_window_along_fork > 0.01))' |\
  mlr --tsv put '$keep_max_300T_sliding_window_along_fork_gt_0pt50 = capitalize(string($max_300T_sliding_window_along_fork > 0.50))' |\
  sed "1i# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>" \
  > "$outputDir"/pauseList_processed

# remove temporary files
rm "$temp_file_1" "$temp_file_2"

# get pause bedgraphs
bash rDNA_pause_coverage.sh "$ref_fasta" "$contig_name" "$num_repeat_units" "$outputDir"/pauseList_processed \
  "$bed_file_with_features" "$outputDir"/pauseBedgraphs \
  "$outputDir"/nascent.bed "$base" "$number" "$pos_boundary" "$outputDir"/parental.bed