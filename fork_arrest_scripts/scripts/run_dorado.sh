#!/bin/bash
#SBATCH --mem=120G
#SBATCH -c 14
#SBATCH -p ei-gpu
#SBATCH --gres=gpu:1
#SBATCH -J baseModCallDorado
#SBATCH --mail-type=END,FAIL
#SBATCH --time=48:00:00

# goal
# ----
# perform reference-free basecalling and modification calling, and align to linear reference genome in one step.
# NOTE: This script is meant to future-proof our pipeline if we switch to a reference-free model.
#       It is not part of the guppy-dnascent v2 pipeline.

# usage
# -----
# sbatch run_dorado.sh $options_json
# options_json: a json file containing the input files and options. Please see the section below for the options.

# json options
# ------------
# The json options file must look like this (fields can be in any order):
# {
#  "pod5_folder": "/path/to/pod5_folder",
#  "output_unfiltered_bam_file": "/path/to/output_unfiltered_bam_file",
#  "output_filtered_bam_file": "/path/to/output_filtered_bam_file",
#  "dorado_sequencing_summary_file": "/path/to/dorado/sequencing_summary.txt",
#  "dorado_dna_model_folder_simplex_hac": "/path/to/dorado_dna_model_folder_simplex_hac",
#  "ref_genome_fasta": "/path/to/reference.fa",
#  "mod_model_file": "/path/to/mod/model/file.pt",
#  "confidential_tmp_dir": "/path/to/confidential/tmp/dir",
#  "only_prim": 1,
#  "min_len": 1000,
#  "min_qual": 10,
#  "adjust_tags_due_to_remora_dorado_ONT_problems": {
#    "fake_mod_bases": "m",
#    "fake_mod_long_names_0": "5mC",
#    "real_mod_bases": "T"
#  },
#  "dorado_version": "0.6.2",
#  "remora_version": "3.1.0"
# }
# Some of these are required and some are optional. Please read the description below. Any parameter with a default
# value is optional. If you don't provide it, the default value will be used.
# If you have any other key-value pairs in the json file, they will be ignored and will not cause an error.

# pod5_folder: folder containing the pod5 files. If you have fast5 files, you can still use this script by passing
#              the fast5 folder as the pod5 folder. We do not endorse this option as the program is much faster
#              with pod5 files and ONT may discontinue fast5 files in the future.
#              So, it is better to convert fast5 to pod5 files and then use this pipeline.
#              See the end of this comment to see how to do this conversion.
#              You cannot specify individual pod5 files or wildcard patterns like /a/b/*.pod5.
#              You must specify a directory with pod5 files and we will look recursively in that directory for pod5
#              files (by recursive we mean we look in all subdirectories as well).
#              Conversion from fast5 to pod5:
#              You can use the pod5 package for this e.g. `pod5 convert fast5 input.fast5 --output output.pod5`.
#              You can load the pod5 package using `source load_package.sh -pod5`.
# output_unfiltered_bam_file: path to the output bam file containing the basecalls and modification calls.
#                             This file will contain reads of all quality and lengths, including primary
#                             and secondary reads.
# output_filtered_bam_file: path to the output bam file containing the basecalls and modification calls after filtering.
#                           Depending on the flags set, this may contain only primary reads, only reads above a
#                           minimum length, and reads above a minimum mapping quality.
# dorado_sequencing_summary_file: path to the output sequencing summary file.
# dorado_dna_model_folder_simplex_hac: path to the dorado model folder that contains the simplex hac model.
#                                      You must have either downloaded this or gotten it from someone.
# ref_genome_fasta: path to the reference genome in fasta format.
# mod_model_file: path to the file produced by remora model training. This will have a name like model_best.pt,
#                 model_final.pt, or model_000020.pt etc. (please don't take these names literally as a future
#                 remora version may name them somewhat differently). You must get this file from someone or generate
#                 it yourself by running our betta pipeline.
# confidential_tmp_dir: path to a temporary directory where intermediate files will be stored. This directory needs
#                       to be confidential as we do not want to share models trained from our betta-remora pipeline.
#                       If this directory does not exist, it will be created. Discuss with someone if you are unsure
#                       about this.
# only_prim: 1 if you want to retain only primary reads, 0 otherwise. Default is 1.
# min_len: minimum length of the read to keep in bp. Default is 1000.
# min_qual: minimum alignment quality of the read to keep. Default is 0.
# adjust_tags_due_to_remora_dorado_ONT_problems: due to small problems in remora and dorado left unfixed/unidentified
#                                                by ONT (but may be fixed in future versions), we have to adjust the
#                                                tags in the bam file before and after modification calling.
#                                                If you don't know what this is, please talk to someone.
#                                                You can omit this field and the associated subfields if you are sure
#                                                that ONT have fixed this problem.
# fake_mod_bases: due to tag problems (see above), we pretend that the modification code is this.
# fake_mod_long_names_0: due to tag problems (see above), we pretend that the modification long name is this.
# real_mod_bases: this is the real modification code. After running dorado with the fake_mod_bases, we will convert
#                 the fake_mod_bases to real_mod_bases in the final bam file. I think the fake long name is not
#                 recorded in the BAM file (except maybe in the header, I don't know), so we don't need to convert
#                 the fake long name back to the correct one.
# dorado_version: version of dorado. Only dorado 0.3.4, 0.6.2, and 0.7.2 are supported. Option is required.
# remora_version: version of remora. Only remora 3.1.0 and 3.2.0 are supported. Option is required.
# NOTE: you've to ensure that the inputs are coordinated. For example, if you trained a modification model with
#       a particular version of remora, dorado, and a particular dorado model file, you must use the same files here.
#       Our pipeline script here will not check for these things!

# outputs
# -------
# A filtered and an unfiltered bam file, both containing the basecalls and modification calls, and a sequencing
# summary file.
# Each BAM file contains a reference-unanchored modification call and an alignment per read.
# These two bits of information put together can give a reference-anchored modification call as well.
# So everything we need is in the BAM files.

# fail on error
set -e

# load jq first
source load_package.sh -jq

# set inputs and output_folder
options_json=${1:-/dev/null}

# ensure that the correct number of arguments are provided
if [ "$#" -ne 1 ];
then
  echo "Illegal number of parameters. Please provide the options json file." >&2
  echo "Usage: sbatch run_dorado.sh \$options_json" >&2
  echo "Please see the script for the required json fields in the options file." >&2
  exit 1
fi

# ensure that the options json file exists
if [ ! -f "$options_json" ];
then
  echo "options json file does not exist: $options_json"
  exit 1
fi

# get options from the json file
pod5_folder=$(jq -r '.pod5_folder // ""' "$options_json")
output_unfiltered_bam_file=$(jq -r '.output_unfiltered_bam_file // ""' "$options_json")
output_filtered_bam_file=$(jq -r '.output_filtered_bam_file // ""' "$options_json")
dorado_sequencing_summary_file=$(jq -r '.dorado_sequencing_summary_file // "/dev/null"' "$options_json")
dorado_dna_model_folder_simplex_hac=$(jq -r '.dorado_dna_model_folder_simplex_hac // ""' "$options_json")
reference=$(jq -r '.ref_genome_fasta // ""' "$options_json")
mod_model_file=$(jq -r '.mod_model_file // ""' "$options_json")
confidential_tmp_dir=$(jq -r '.confidential_tmp_dir // ""' "$options_json")
onlyPrim=$(jq -r '.only_prim // 1' "$options_json")
minLen=$(jq -r '.min_len // 1000' "$options_json")
minQual=$(jq -r '.min_qual // 0' "$options_json")
fake_mod_bases=$(jq -r '.adjust_tags_due_to_remora_dorado_ONT_problems.fake_mod_bases // ""' "$options_json")
fake_mod_long_names_0=$(jq -r '.adjust_tags_due_to_remora_dorado_ONT_problems.fake_mod_long_names_0 // ""'\
 "$options_json")
real_mod_bases=$(jq -r '.adjust_tags_due_to_remora_dorado_ONT_problems.real_mod_bases // ""' "$options_json")
dorado_version=$(jq -r '.dorado_version // ""' "$options_json")
remora_version=$(jq -r '.remora_version // ""' "$options_json")

# set adjustTags to 1 if fake_mod_bases, fake_mod_long_names_0, and real_mod_bases are not empty
# or set it to 0
if [ "$fake_mod_bases" != "" ] && [ "$fake_mod_long_names_0" != "" ] && [ "$real_mod_bases" != "" ];
then
  adjustTags=1
else
  adjustTags=0
fi

# load samtools, modkit, and python
source load_package.sh -samtools -modkit -python

# load remora
if [ "$remora_version" == "3.1.0" ];
then
  source load_package.sh -remora_3_1_0
elif [ "$remora_version" == "3.2.0" ];
then
  source load_package.sh -remora_3_2_0
else
  echo "remora version $remora_version is not supported. Only remora 3.1.0 and 3.2.0 are supported." >&2
  exit 1
fi

# load dorado
if [ "$dorado_version" == "0.6.2" ];
then
  source load_package.sh -dorado_0_6_2
elif [ "$dorado_version" == "0.7.2" ];
then
  source load_package.sh -dorado_0_7_2
elif [ "$dorado_version" == "0.3.4" ];
then
  source load_package.sh -dorado_0_3_4
else
  echo "dorado version $dorado_version is not supported. Only dorado 0.6.2 and 0.7.2 are supported." >&2
  exit 1
fi

# ensure that the input folder exists
if [ "$pod5_folder" == "" ] || [ ! -d "$pod5_folder" ];
then
  echo "pod5 folder does not exist: $pod5_folder" >&2
  exit 1
fi

# complain if any of the required options are missing
if [ ! -d "$dorado_dna_model_folder_simplex_hac" ];
then
  echo "dorado dna model folder does not exist: $dorado_dna_model_folder_simplex_hac" >&2
  exit 1
fi

if [ ! -f "$reference" ];
then
  echo "reference genome does not exist: $reference"  >&2
  exit 1
fi

if [ ! -f "$mod_model_file" ];
then
  echo "modification model file does not exist: $mod_model_file"  >&2
  exit 1
fi

# set temporary directory. throw an error if confidential_tmp_dir is empty
if [ "$confidential_tmp_dir" == "" ];
then
  echo "confidential_tmp_dir is empty. Please provide a confidential temporary directory." >&2
  exit 1
fi
tmpDir="$confidential_tmp_dir"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# export the remora model
model_folder="$tmpDir"/mod_model_folder
mkdir -p "$model_folder"
remora model export "$mod_model_file" "$model_folder"

# adjust the tag of the config.toml file
if [ "$adjustTags" -eq 1 ];
then
  python adjust_config_toml.py "$model_folder"/config.toml "$fake_mod_bases" "$fake_mod_long_names_0"
fi

# make any necessary directories
mkdir -p "$(dirname "$output_unfiltered_bam_file")"
mkdir -p "$(dirname "$output_filtered_bam_file")"
mkdir -p "$(dirname "$dorado_sequencing_summary_file")"

# run dorado
dorado  \
      basecaller \
      "$dorado_dna_model_folder_simplex_hac" \
      "$pod5_folder" \
      --recursive \
      --reference "$reference" \
      --modified-bases-models  "$model_folder" \
      --emit-moves \
      --modified-bases-threshold 0 \
     > "$tmpDir"/mod_calls.bam

if [ "$adjustTags" -eq 1 ];
then
  # adjust the tags in the bam file
  modkit adjust-mods --threads 10 --convert "$fake_mod_bases" "$real_mod_bases" \
    "$tmpDir"/mod_calls.bam "$tmpDir"/mod_calls.adjusted.bam
else
  # just rename the file
   mv "$tmpDir"/mod_calls.bam "$tmpDir"/mod_calls.adjusted.bam
fi

# sort and index the reads
samtools sort -@ 16 -o "$output_unfiltered_bam_file" \
  -T "$tmpDir" "$tmpDir"/mod_calls.adjusted.bam
samtools index "$output_unfiltered_bam_file"

# get a sequencing summary file
dorado summary "$output_unfiltered_bam_file" > "$dorado_sequencing_summary_file"

# run samtools
if [ "$onlyPrim" -eq 1 ];
then
  # get flags ready for primary read filtering so only primary reads are allowed through
  secFlag=256
  supplFlag=2048
  invPrimFlag=$((secFlag + supplFlag))
else
  # no flag-based filtering
  invPrimFlag=00
fi

# filter the reads
samtools view -e "rlen>=$minLen&&mapq>=$minQual" -b -F $invPrimFlag \
  -o "$tmpDir"/mod_calls_filtered.bam "$output_unfiltered_bam_file";

# sort and index the reads
samtools sort -@ 16 -o "$output_filtered_bam_file" \
  -T "$tmpDir" "$tmpDir"/mod_calls_filtered.bam
samtools index "$output_filtered_bam_file"

# remove the temporary folder
rm -rf "$tmpDir"
