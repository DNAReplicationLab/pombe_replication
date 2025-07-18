#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J pauseAnalysis
#SBATCH --mail-type=END,FAIL
#SBATCH --time=64:59:59
#SBATCH --constraint=""

# Goal
# -----
# Given pauses at genome locations in a nanopore dataset of budding yeast with sigmoidal BrdU fall off, calculate
# enrichment of pauses at different features and compare them with a sophisticated calculated null hypothesis.

# Prerequisites
# -------------
# Before running this script, you must have run the pause detection pipeline (get_sgm_dnascent_pauses.sh)
# on a dataset and you must have locations of bed files and input options corresponding to different features
# like tRNAs, genes etc. on the reference genome. For more details, see the section usage below.

# usage
#------
# bash pipeline_sgm_pauses_winnow_and_analyse.sh <dataset_json_options> <calculation_json_options>
# dataset_json_options: json file with options for loading dataset, see the section below for details
# calculation_json_options: json file with options for calculating pause enrichment, see the section below for details

# contents of json input files
# ----------------------------
# Introduction:
# A json file is made up of a bunch of key-value pairs. The key is a string, and the value can be a string, number,
# array, or another json object. The key-value pairs are separated by commas, and the entire json object is surrounded
# by curly braces. For example, the following is a valid json object and can be saved in a file called "example.json":
# {
#   "key1": "value1",
#   "key2": 2,
#   "key3": [1, 2, 3],
#   "key4": {
#     "key5": "value5"
#   }
# }
# NOTE: We have set it up so that we can do variable substitution in the json files.
#       Consider the following two key-value pairs:
#       "key1": "value1", "key2": "value2_${key1}"
#       The value of key2 will be "value2_value1" after variable substitution.
#       If you don't understand what this means, don't worry about it and don't use it.

# Details of dataset_json_options
# -------------------------------
# dataset_json_options should contain the following keys and appropriate values:
#  mod_bam: path to the modified bam file with analogue probabilities vs position per read
#  forksense_dir: directory with fork sense files produced by dnascent containing left fork, right fork, origin,
#                 termination information. NOTE: If you've called pauses on forks after filtering them on some
#                 constraint e.g. fork length > 5 kb, alignment length > 30 kb or something like that, then
#                 you must use a directory with fork sense files with these filters applied.
#                 The files must be named with the standard forkSense names i.e. leftForks_DNAscent_forkSense.bed,
#                 rightForks_DNAscent_forkSense.bed etc. So if you have files named differently, you must rename them.
#  input_pause_file: pause file called by our pause-detection pipeline, in tab-separated value format
#  dataset: a string which labels the dataset from which the input files were derived e.g. "20180926_JK_ONT_SC_wt40_a07aba2"
#  seq_summary_file: path to the sequencing summary file for this dataset, produced by the basecaller.
#  fasta_file: path to the fasta file containing the genome sequence. must have a corresponding .fai file in the
#              same directory. You can use sacCer3.fa here, no need to worry about chrM or that there are only
#              two copies of rDNA in the reference genome.
#  output_dir: path to the directory where output files will be written. This directory will be created if it doesn't
#              exist. Choose an empty/new directory unless you are sure there will be no overwriting of files.

# Example:
# {
#    "mod_bam": "/path/to/analogue.mod.bam",
#    "forksense_dir": "/path/to/forkSenseOverallBedgraphs",
#    "dataset": "20230824_JD_ONT_SC_JD01744uM_abe3689",
#    "input_pause_file": "/path/to/pause",
#    "seq_summary_file": "/path/to/sequencing_summary.txt",
#    "fasta_file": "/path/to/sacCer3.fa",
#    "output_dir": "/path/to/output/directory"
# }

# Details of calculation_json_options
# -----------------------------------
# calculation_json_options should contain the following keys and appropriate values:
# NOTE: the pipeline_go_* keys are optional. If you want to run a particular stage of the pipeline, then set the
#       corresponding key-value pair to 1. If you don't want to run a particular stage of the pipeline, then
#       you can omit the corresponding key-value pair or set it to 0.
# pipeline_go_calculate_sensitivity: switch on the pipeline stage to calculate pause sensitivity
# pipeline_go_plot_sensitivity_vs_trep: switch on the pipeline stage to plot pause sensitivity vs median replication
#                                       timing. uses the t_rep_file described below
# pipeline_go_calculate_enrichment: switch on the pipeline stage to calculate pause enrichment, uses the feature_files
#                                   described below
# do_not_redo_sensitivity_calculations: (default 0) set to 1 if you do not want to overwrite sensitivity files.
#                                       i.e. this option runs the sensitivity calculation only if the associated
#                                       output files do not exist.
# analysis_label: a short string which labels the analysis done here e.g.: a date like "28sep23"
# log_file: (default /dev/null i.e. no logging) path to the log file where the steps of the pipeline will be logged.
# t_rep_file: path to the t_rep file containing the Median replication timing information in 1 kb bins.
#             This file is not produced by us, but is from a different publication.
#             You should get this from someone. A few lines of this file look like
#             (space-separated, don't assume the below is real data):
# chrI 3000 3999 50.92
# chrI 3999 4999 51.72
# chrI 4999 6000 53.28
# feature_files: - A list of json files, each of which contains the keys and appropriate values in the list below.
#           - For more details, see the script that receives these json files as inputs in our code below.
#           - We are brief here as we don't want to repeat comments from another script.
#           - Many more input options are possible; please refer to the script that receives these json files as inputs.
#             e.g.: no_jobs_for_head_on_and_co_directional
#           - The idea is that we split a feature e.g. genes into equal bins by length and into equal bins
#             by a numeric measurement like transcription rate.
#           - Note that some of these options are optional e.g.: if you only want to split by length and don't want
#             to split by both length and some numeric value, then you can drop some options
#   feature: a string which labels the feature e.g. "tRNA"
#   input_bed_file: path to the bed file containing the locations of the feature.
#   features: instead of using feature and input_bed_file, you can use the features list to specify multiple features.
#             Each feature in the list should have the keys feature and input_bed_file.
#             NOTE: You can use either this or feature and input_bed_file, not both.
#             e.g.: "features": [ { "feature": "TSS", "input_bed_file": "/path/to/TSS.bed" },
#                                 { "feature": "TES", "input_bed_file": "/path/to/TES.bed" } ]
#             If you use this option, you can set the feature to "auto" to get the feature name from the bed file.
#             The auto name option will work only for some datasets, so it is recommended to use the feature name
#             explicitly.
#   column_of_interest: a number which labels a numeric column of the bed file e.g.: transcription rate.
#   num_split_by_value: a number which labels the number of bins to split the column of interest into.
#   split_by_length_bp: a number which labels the length of each bin in base pairs when the feature is split by length.
#   align_by: e.g. "head-to-tail", the way in which features are aligned.
#   genome_size_bp_optional: size of the genome in base pairs, must include rDNA with appropriate copy number,
#                            must exclude chrM.

# Example
# {
#   "t_rep_file": "/path/to/t_rep.bedgraph",
#   "feature_files": [ "/path/to/feature1.json", "/path/to/feature2.json" ],
#    "log_file": "/path/to/log/file"
# }

# outputs
# -------
# Many output files are sent to the output directory. They fall in a few categories:
# pause sensitivity bedgraphs: These contain calculations corresponding to our null hypothesis i.e.
#                              to establish if pauses are enriched/depleted more than expected at random at a
#                              given a genomic region, we need to have a null hypothesis of how many pauses to expect
#                              at that region. These bedgraph files contain a signal vs genomic position.
#                              Summing this signal over a region gives the expected number of pauses at that region.
#                              For more details, refer to the corresponding script below.
# feature directories: We calculate enrichment and depletion at different genomic features specified in the input
#                      calculation_json_options file. For each feature, we create a directory with the name of the
#                      feature in the output directory. Each directory can contain 100s or 1000s of files.
#                      For more information, refer to the corresponding script below.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -jq -python

# load configuration
source config.sh

# print node information
bash print_node_information.sh

# create the function for extracting PIPELINE_GO variables
get_pipeline_go(){

  PIPELINE_GO=$(jq -r '.'"$2"'' "$1")
  if [ "$PIPELINE_GO" != "0" ] && [ "$PIPELINE_GO" != "1" ]; then
      echo "pipeline_go variables should be 0 or 1, so setting it to 0" >&2
      echo 0;
  else
      echo "$PIPELINE_GO";
  fi

}

# assign arguments to variables
dataset_json_options=${1:-}
calculation_json_options=${2:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 <dataset_json_options> <calculation_json_options>"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# check that the two json files exist
if [ ! -f "$dataset_json_options" ]; then
    >&2 echo "ERROR: dataset_json_options file does not exist"
    exit 1;
fi

if [ ! -f "$calculation_json_options" ]; then
    >&2 echo "ERROR: calculation_json_options file does not exist"
    exit 1;
fi

# make a temporary directory
tmpDir=$(mktemp -p "${config[scratchDir]:-}"/tmp -d)
mkdir -p "$tmpDir"

# expand parameters in input json files if need be
< "$dataset_json_options" python expand_variables_in_json.py > "$tmpDir"/dataset.json
< "$calculation_json_options" python expand_variables_in_json.py > "$tmpDir"/calculation.json

# make output directory if it doesn't exist
output_dir=$(jq -r '.output_dir // "/tmp"' "$tmpDir"/dataset.json)
mkdir -p "$output_dir"
output_dir=$(realpath "$output_dir")

# load files from json inputs
pause_file=$(jq -r '.input_pause_file' "$tmpDir"/dataset.json)
fasta_file=$(jq -r '.fasta_file' "$tmpDir"/dataset.json)
t_rep=$(jq -r '.t_rep_file' "$tmpDir"/calculation.json)
seq_summary_file=$(jq -r '.seq_summary_file' "$tmpDir"/dataset.json)
analysis_label=$(jq -r '.analysis_label' "$tmpDir"/calculation.json)
log_file=$(jq -r '.log_file // "/dev/null"' "$tmpDir"/calculation.json)

# load options from json input
is_do_not_redo_sens_calc=$(jq -r '.do_not_redo_sensitivity_calculations // 0' "$tmpDir"/calculation.json)

# Step 0
# ======
# begin logging after making the required directory for the log file
mkdir -p "$(realpath "$(dirname "$log_file")")"
{
  echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Beginning pause analysis calculations"
  echo "If slurm job, the job id is $SLURM_JOB_ID"
  echo "Finished parsing inputs"
} >> "$log_file"

# Step 1: Calculate pause sensitivities for pauses restricted to those with 99% CI in pause site less than +-100 bases
# ====================================================================================================================

# check that input pause file exists
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: input_pause_file $pause_file does not exist"
    exit 1;
fi

PIPELINE_GO=$(get_pipeline_go "$tmpDir"/calculation.json "pipeline_go_calculate_sensitivity")
if [ "$PIPELINE_GO" -eq 1 ];
then
  # check that fasta file exists
  if [ ! -f "$fasta_file" ]; then
      >&2 echo "ERROR: fasta_file file does not exist"
      exit 1;
  fi

  # set the -w flag if the user does not want to overwrite sensitivity files
  sens_flag=""
  if [ "$is_do_not_redo_sens_calc" -eq 1 ]; then
    sens_flag="-w"
  fi

  echo "Submitting a pause sensitivity job for the restricted pause file $pause_file" >> "$log_file"
  echo "NOTE: we run sbatch in wait mode here, so logs will continue after the job finishes" >> "$log_file"
  pause_sens_job=$(
    sbatch --wait run_get_sgm_pause_sensitivities.sh $sens_flag "$pause_file"\
      keep_fp,keep_fl,keep_fs\
      "$fasta_file"\
      "$pause_file"\
      10 ); # window size 10, decimate_T: default
  echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Pause sensitivity job submitted: $pause_sens_job" >> "$log_file"
fi

echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Will use pause file $pause_file for subsequent steps" >> "$log_file"

# Step 2: Plot median replication time versus pause sensitivity for the newly calculated (restricted) pause sensitivity
# =====================================================================================================================
PIPELINE_GO=$(get_pipeline_go "$tmpDir"/calculation.json "pipeline_go_plot_sensitivity_vs_trep")
if [ "$PIPELINE_GO" -eq 1 ];
then

  # check if output directory is set to /tmp, then exit
  if [ "$output_dir" = "/tmp" ]; then
    >&2 echo "ERROR: output_dir is set to the default /tmp, which is not allowed in this step"
    exit 1;
  fi

  pause_sens="$pause_file".pause_sensitivity.all.bedgraph
  op_prefix="$output_dir"/pause_sens

  # check that the required input files exist
  if [ ! -f "$t_rep" ]; then
    >&2 echo "ERROR: t_rep file does not exist"
    exit 1;
  fi

  if [ ! -f "$fasta_file".fai ]; then
    >&2 echo "ERROR: fasta_file.fai file does not exist"
    exit 1;
  fi

  if [ ! -f "$seq_summary_file" ]; then
    >&2 echo "ERROR: seq_summary_file file does not exist"
    exit 1;
  fi

  plot_trep_vs_pause_sens_job=$(
    sbatch -p ei-medium --time=02:00:00 calculate_and_plot_trep_vs_pause_sens.sh\
      "$t_rep" "$pause_sens" "$fasta_file".fai "$seq_summary_file" "$op_prefix"
  );
  echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Pause sensitivity vs trep plotting job submitted: $plot_trep_vs_pause_sens_job" \
    >> "$log_file"
else
  echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Will not calculate pause sensitivity vs trep and plot it" >> "$log_file"
fi

# Step 3: Perform pause enrichment calculation per feature
# =========================================================
PIPELINE_GO=$(get_pipeline_go "$tmpDir"/calculation.json "pipeline_go_calculate_enrichment")
if [ "$PIPELINE_GO" -eq 1 ];
then

  # check if output directory is set to /tmp, then exit
  if [ "$output_dir" = "/tmp" ]; then
    >&2 echo "ERROR: output_dir is set to the default /tmp, which is not allowed in this step"
    exit 1;
  fi

  # now, perform enrichment calculation per feature
  n_feature_files=$(jq -r '.feature_files | length' "$tmpDir"/calculation.json)

  for count in $(seq 0 "$((n_feature_files-1))"); do

    # first, extract the feature file contents to a temporary file
    tmp_feature_file=$(mktemp --tmpdir="$tmpDir" feature_"$count"_XXX.json)
    cp "$(jq -r '.feature_files['"$count"']' "$tmpDir"/calculation.json)" "$tmp_feature_file"

    # check if this file contains the keys of "features" or ("feature" and "input_bed_file").
    # If it contains "features", then we need to expand the list of features into separate files.
    # If it contains "feature" and "input_bed_file", then we don't need to expand the list of features.
    # If it contains both, then we throw an error.
    is_appropriate=$(jq '
      def xor($a; $b):
        ($a or $b) and (($a and $b) | not)
      ;
      xor(has("features"); has("feature") and has("input_bed_file"))
      ' "$tmp_feature_file")

    if [ "$is_appropriate" = "false" ]; then
      echo "ERROR: feature file should contain either 'features' or ('feature' and 'input_bed_file')" >&2
      exit 1;
    elif [ "$(< "$tmp_feature_file" jq -r 'has("features")')" = "true" ]; then

        # expand the features into separate files
        len_features=$(< "$tmp_feature_file" jq -r '.features | length')
        for count_2 in $(seq 1 "$len_features"); do

          feature_file=$(mktemp --tmpdir="$tmpDir" feature_"$count"_XXX.json)
          current_feature=$(< "$tmp_feature_file" jq -r '.features['"$count_2-1"'].feature // "NA"')
          current_input_bed_file=$(< "$tmp_feature_file" jq -r '.features['"$count_2-1"'].input_bed_file // "NA"')

          # check that the current feature is not empty
          if [ "$current_feature" = "NA" ]; then
            >&2 echo "ERROR: feature is empty"
            exit 1;
          fi

          # if feature is set to auto, then get the feature from the bed file
          if [ "$current_feature" = "auto" ]; then
            current_feature=$(python get_feature_name.py "$(basename "$current_input_bed_file")")
          fi

          # create a new feature file retaining the other keys and values but updating the feature and input_bed_file
          jq --arg feature "$current_feature" --arg input_bed_file "$current_input_bed_file" \
            '.feature = $feature | .input_bed_file = $input_bed_file | del(.features)' "$tmp_feature_file" \
              > "$feature_file"

        done

        rm "$tmp_feature_file"
    fi

    for feature_file in "$tmpDir"/feature_"$count"*; do

      # get a random string to coordinate jobs
      rand_string=$(openssl rand -hex 2)

      # create a temporary calculation options file
      calculation_json_options="$tmpDir"/calculation_json_options_"$rand_string".json

      # first, combine the feature file and our dataset options
      # NOTE: If any fields are repeated, the latest one is retained and others are discarded
      {
        echo "{"
        jq . "$feature_file" | sed '1d' | sed '$d' | sed '$s/$/,/'
        jq . "$tmpDir"/dataset.json | sed '1d' | sed '$d' | sed '$s/$/,/'
        echo "\"pause_file\": \"$pause_file\","
        echo "\"pause_sensitivity_bedgraphs_prefix_optional\": \"${pause_file}.pause_sensitivity\","
        echo "\"input_pause_file\": \"NA\","
        echo "\"analysis_label\": \"$analysis_label\","
        echo "\"random_string_job_name\": \"bed_pause_$rand_string\","
        echo "\"log_file\": \"$log_file\","
        echo "\"output_dir\": \"$output_dir/\${feature}\","
        echo "\"op_dir\": \"\${output_dir}\","
        echo "\"output_prefix\": \"\${output_dir}/bed_files/bed\""
        echo "}"
      } | jq . | python expand_variables_in_json.py > "$calculation_json_options"

      # split bed file and run many bed region pause report
      launchJob=$(sbatch -p ei-medium --time=06:00:00 --wait -c 1 --mem=10G \
        split_bed_and_run_many_bed_region_pause_report.sh "$calculation_json_options");
      jid=${launchJob##* }
      echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Split bed and generate pause report job with ID $jid submitted and finished."\
        >> "$log_file"

      # plot results
      pause_report_directory=$(jq -r '.output_dir' "$calculation_json_options")
      plot_output_directory="$pause_report_directory"/bed_split_pause_summary
      launchJob=$(sbatch --dependency=singleton --job-name=bed_pause_"$rand_string" \
        run_plot_pause_profiles_from_pause_reports_from_split_bed.sh "$pause_report_directory" "$plot_output_directory");
      jid=${launchJob##* }
      echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Plot pause count results job with ID $jid submitted." >> "$log_file"

    done

  done
else
  echo "$(date +\[%Y_%m_%d_%H:%M:%S\]) Will not calculate feature enrichment" >> "$log_file"
fi