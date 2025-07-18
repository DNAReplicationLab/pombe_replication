#!/bin/bash

# goal
# -----
# Go from nanopore data (this could be raw current data or basecalled data) to modification calls and processing
# i.e. perform modification detection and some preliminary qualitative and quantitative analysis.
# There are three modes here: reference-anchored R9 dnascent 2, reference-anchored R10 dnascent 4,
# and reference-free R10 dorado. The mode is chosen based on the parameters in the input json file.

# overall notes
# -------------
# * Although the R10-dorado-ref-free mode is called "ref-free", the adjective refers only to the modification
#   calling process. We do an alignment to the reference genome that is independent of the modification calling.
#   Using the modification calls and the alignment information one can infer a reference-anchored modification profile
#   as well. This is why we require a reference genome for this mode, and any tool that needs
#   reference-anchored basecalling can accept these bam files as input.
# * Although we say this pipeline can run reference-anchored R10 dnascent 4, the pipeline is only
#   equipped to deal with one analogue (this may change in future) whereas dnascent 4 calls both BrdU
#   and EdU probabilities per reference position per read. So, we just discard the EdU calls and
#   convert to a .detect file in the DNAscent version 2 format. The detect file in the original format
#   output by dnascent 4 is stored in a suitably named backup file if the EdU calls are needed in the future.
# * DNAscent mode means DNAscent detect v4 or v2 will be used to call modifications.
#   DNAscent also includes a package called forkSense to call forks using BrdU gradients.
#   This package is downstream of modification calling.
#   This package may be used in any of the three modes, so please do not be confused if you see
#   DNAscent forkSense being used in the R10-dorado-ref-free mode.
# * DNAscent detect may have trouble if you convert pod5 from the sequencer -> fast5 but use the basecalled fastq files
#   generated during sequencing. If you are running into this problem, just re-basecall the fast5 files using this
#   pipeline and you should be fine.

# usage
# -----
# bash pipeline_fast5_pod5_to_dnascent_dorado.sh $input_json
# - no need to use sbatch as this script just submits other jobs and does not do any other computation.
#   So the script will be fast enough in just plain bash.
# - input_json is the path to the json file with fields required by the scripts called here.
# - the contents of input_json look like the following, with sample inputs (fields can be in any order):
# {
# "ref_genome_fasta": "/path/to/ref_genome.fasta",
# "min_len": 1000,
# "min_qual": 20,
# "only_prim": 1,
# "dataset_prefix": "20190205_ST_ONT_SC_wt",
# "force_no_commit_str": 1,
# "fast5_dir": "/path/to/fast5_dir",
# "fastq_dir": "/path/to/fastq_dir",
# "fastq_file": "/path/to/fastq_file",
# "guppy_config_file": "dna_r9.4.1_450bps_hac.cfg",
# "guppy_sequencing_summary_file": "/path/to/guppy/sequencing_summary.txt",
# "plot_bed_file": "/path/to/plot_bed_file",
# "hdf5_plugin_path": "/path/to/hdf5/plugin",
# "make_db_new_fast5_format": 1,
# "use_dnascent_version_4": 1,
# "guppy_version": "5.0.7",
# "use_R10_dorado_ref_free_model": 1,
# "pod5_folder": "/path/to/pod5/files",
# "dorado_dna_model_folder_simplex_hac": "/path/to/dorado/model/folder",
# "mod_model_file": "/path/to/mod/model/file.pt",
# "confidential_tmp_dir": "/path/to/confidential/tmp/dir",
# "dorado_sequencing_summary_file": "/path/to/dorado/sequencing_summary.txt",
# "adjust_tags_due_to_remora_dorado_ONT_problems": {
#    "fake_mod_bases": "m",
#    "fake_mod_long_names_0": "5mC",
#    "real_mod_bases": "T"
# },
# "dorado_version": "0.6.2",
# "remora_version": "3.1.0",
# "pipeline_go_R10_ref_free_basecall_modcall_align": 1,
# "pipeline_go_R10_ref_free_convert_to_dnascent2_format": 1,
# "pipeline_go_basecall": 1,
# "pipeline_go_align": 1,
# "pipeline_go_filter_sort_before_dnascent": 1,
# "pipeline_go_pycoQC": 1,
# "pipeline_go_make_db": 1,
# "pipeline_go_dnascent_index": 1,
# "pipeline_go_detect": 1,
# "pipeline_go_make_modbam": 1,
# "pipeline_go_annotate_modbam": 1,
# "pipeline_go_if_dnascent4_then_get_dnascent2_format": 1,
# "pipeline_go_forkSense": 1,
# "pipeline_go_dnascent_qc": 1,
# "pipeline_go_make_forkSense_modbam": 1,
# "pipeline_go_plot_random_nascent_reads": 1,
# "pipeline_go_plot_random_nascent_reads_forkSense": 1,
# "pipeline_go_plot_random_nascent_reads_collate_pdf": 1,
# "pipeline_go_plot_filtered_random_nascent_reads": 1
# }
# To know what the fields mean, please look at the corresponding section in the script.
# We provide a brief description here:
# MODES: * the R10-dnascent4-ref-anchored mode is chosen if use_dnascent_version_4 is set to 1.
#        * the R10-dorado-ref-free mode is chosen if use_R10_dorado_ref_free_model is set to 1.
#        * If the above two are set to 0 or not set, then the R9-dnascent2-ref-anchored mode is chosen.
#        * The required inputs for each mode are different, so please read the corresponding sections in the script.
#          and set them appropriately. The script will complain if you don't set them correctly.
#        * Parameters descriptions below specify the mode they are needed for e.g. "R10-dnascent4-ref-anchored",
#          "R10-dorado-ref-free", "see run_dorado.sh"/"used in dorado" (the implication is that it is needed
#          for R10-dorado-ref-free). If a mode is not mentioned, then the parameter is needed for all modes.
#          The script will check if the parameters are set correctly for the mode chosen but it is possible
#          that the script may not catch all errors, so please set them carefully after reading the script.
# CAUTION: Do not skip reading all the comments in the script before setting variables.
#          There are some dependencies between variables or other details that you should be aware of.
#          For e.g. if you do not set pipeline_go_basecall to 1, then the guppy_config_file option is irrelevant.
# * ref_genome_fasta: path to the reference genome fasta file.
# * min_len: minimum length of reads to be used for modification calling (in bp).
# * min_qual: minimum quality of reads to be used for modification calling.
# * only_prim: whether to use only primary reads (1) or not (0) in modification calling.
# * dataset_prefix: prefix of the dataset. This is used to name the output files.
#                   see the corresponding section in the script to learn our naming convention.
# * force_no_commit_str: (default 0) usually _commitstr is appended to the dataset_prefix to create a unique name
#                        for the output. If you want to force the output to be named dataset_prefix without the suffix,
#                        then set this to 1.
# * fast5_dir: path to the directory containing fast5 files, omit if using R10-dorado-ref-free.
# * fastq_dir: path to the directory containing fastq files, omit if using R10-dorado-ref-free.
#              leave this blank ("") or unspecified if you want to generate these files in the pipeline here.
#              for more info, see the corresponding section in the script.
# * fastq_file: path to the fastq file(s), omit if using R10-dorado-ref-free.
#               leave this blank ("") or unspecified if you want to generate these files in the pipeline here.
#               for more info, see the corresponding section in the script.
# * guppy_config_file: the guppy config file without the path e.g. dna_r9.4.1_450bps_hac.cfg, omit if
#                      R10-dorado-ref-free.
#                      This is only needed if pipeline_go_basecall is set to 1, otherwise it is irrelevant but it
#                      is a good practice to set it for record keeping. This is the model used by guppy for basecalling.
#                      The model can be changed according to which pore chemistry is used, what accuracy is needed etc.
# * guppy_sequencing_summary_file: path to the guppy sequencing summary file, this is already available or will be
#                                  created during guppy basecalling. leave this blank "" or unspecified if you want
#                                  to generate these files in the pipeline here. Omit if using R10-dorado-ref-free.
#                                  NOTE: a sequencing summary file is normally generated by the nanopore device when
#                                  on-board basecalling is done. If you are re-basecalling using the pipeline here,
#                                  then you want to create a new sequencing summary file associated with this
#                                  basecalling. In this scenario, you do not want to overwrite the sequencing summary
#                                  generated by on-board basecalling. So, you should set this parameter to a new file,
#                                  or leave it blank so that a new file with some default name is created by us in the
#                                  output directory.
# * dorado_sequencing_summary_file: same as the guppy sequencing summary file, but for dorado. Same instructions
#                                   as above but for dorado basecalling. Only used if using R10-dorado-ref-free.
#                                  NOTE: a sequencing summary file is normally generated by the nanopore device when
#                                  on-board basecalling is done. If you are re-basecalling using the pipeline here,
#                                  then you want to create a new sequencing summary file associated with this
#                                  basecalling. In this scenario, you do not want to overwrite the sequencing summary
#                                  generated by on-board basecalling. So, you should set this parameter to a new file,
#                                  or leave it blank so that a new file with some default name is created by us in the
#                                  output directory.
# * plot_bed_file: (can be omitted) path to the bed file containing annotations for the mod bam e.g. gene locations,
#                  tRNA locations etc. Whatever the user desires.
#                  For more details, see the comments in annotate_mod_bam_with_bed.sh.
# * hdf5_plugin_path: path to the hdf5 plugin needed by dnascent. Omit if using R10-dorado-ref-free.
#                     If the plugin is /a/b/c/libvbz_hdf_plugin.so then set this to /a/b/c.
# * make_db_new_fast5_format: (default 0) set this to 1 to tell the db creation script that the fast5 files are in
#                             the new format. Set this to zero or skip if the fast5 files are in the old format.
#                             Irrelevant if pipeline_go_make_db is set to 0. Omit if using R10-dorado-ref-free.
# * use_dnascent_version_4: (default 0) set this to 1 if you want to use dnascent version 4. Set this to zero or skip
#                           if you want to use version 2. Irrelevant if all three of pipeline_go_detect,
#                           pipeline_go_if_dnascent4_then_get_dnascent2_format, and
#                           pipeline_go_dnascent_index are set to 0. If you set this to 1 and the basecalling step is
#                           set to 1, then you must set the guppy config file to an r10 model and the guppy_version
#                           to 6.5.7 (If you don't know what this means or how to do this; please talk to someone
#                           or google it). If you set use_dnascent_version_4 to 1 and you set pipeline_go_detect to 1,
#                           then you must set pipeline_go_if_dnascent4_then_get_dnascent2_format to 1.
#                           Omit if using R10-dorado-ref-free or R9-dnascent2-ref-anchored.
# * guppy_version: (default 5.0.7) the version of guppy used for basecalling. Set this to 6.5.7 if you are using
#                  dnascent version 4. The default is 5.0.7 for historical reasons. For now, you cannot use
#                  any other version of guppy, but this may change as the pipeline script is updated.
#                  Omit if using R10-dorado-ref-free.
# * use_R10_dorado_ref_free_model: (default 0) set this to 1 if you want to use the R10 dorado ref-free model.
# * pod5_folder: path to the directory containing pod5 files. Please see run_dorado.sh for more details, including
#                what to do if you have fast5 files from R10 but not pod5.
# * dorado_dna_model_folder_simplex_hac: path to the dorado model folder that contains the simplex hac model.
#                                        You must have either downloaded this or gotten it from someone.
# * mod_model_file: path to the remora modification model file. See run_dorado.sh for more details.
# * confidential_tmp_dir: path to a directory where temporary files will be stored, only used if R10-dorado-ref-free.
#                         This directory should be confidential because we do not want our modification models
#                         generated by the betta-remora pipeline to be shared at this moment.
# * adjust_tags_due_to_remora_dorado_ONT_problems: We need to set some fake tags due to (small) problems with ONT
#                                                  software. Please see run_dorado.sh for more details.
#                                                  This is a json object with a few fields.
# * dorado_version: version of dorado. See run_dorado.sh for more details and for what values are allowed.
# * remora_version: version of remora. See run_dorado.sh for more details and for what values are allowed.
# * pipeline_go options: set them to 1 if you want to run the corresponding step. any omitted option defaults to 0.
#                        To know what each stage does, please read the code corresponding to that stage in the script.

# initial steps to prepare for execution
# ======================================
set -e

# load git labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -jq

# load some config information
source config.sh

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

# check inputs
# ============

# check that there is one command line argument
if [ "$#" -ne 1 ]; then
    echo "usage: bash pipeline_fast5_to_dnascent.sh input.json" >&2
    echo "see the script for details on what should be in the input.json file" >&2
    exit 1
fi

# check that the input json file exists
if [ ! -f "$1" ]; then
    echo "input json file does not exist" >&2
    exit 1
fi

input_json="$1"

# declare all directories, filenames, some parameters
# ===================================================

# NOTE: none of these paths must have spaces in them.

# fasta reference genome
# ----------------------
refGenomeFasta=$(jq -r '.ref_genome_fasta' "$input_json")

# check that the reference genome exists
if [ ! -f "$refGenomeFasta" ]; then
    echo "reference genome $refGenomeFasta does not exist"
    exit 1
fi

# some read/alignment filtering parameters
# -----------------------------------------
# after minimap is run, we run another samtools command
# to filter reads. these parameters are fed to that command,
# and to dnascent.
minLen=$(jq -r '.min_len' "$input_json")         # lower threshold of read length in bases
minQual=$(jq -r '.min_qual' "$input_json")       # lower threshold of read mapping quality
onlyPrim=$(jq -r '.only_prim' "$input_json")     # set this flag to one if only primary reads should be retained
               # warning: if flag is not set, there may be multiple alignments per read id.
               #          so unset only if you know what you are doing.

# output directory
# ----------------
# NOTE: if analyzing dataset made within the CAN group,
#       the name of output directory must have the format
#       str1_str2_str3_str4_str5_str6
#       str1=date data was generated in yyyymmdd format
#       str2=experimenter's initials in two or three letters
#       str3=sequencing method, probably ONT for (Oxford) nanopore
#       str4=two letters of the species
#       str5=some human readable text WITHOUT special characters like underscores or spaces
#       str6=automatically generated by the system, the first few characters of the git commit hash
# NOTE: if analyzing an external dataset, use a suitable name.
#       if analyzing a dataset associated with a paper, use str1_str2_str3_str4
#       where str1=first author's surname
#             str2=short journal name, like pnas or ncomms
#             str3=year in yy format
#             str4=use commit string generated by the system, available in the COMMITSTR variable
# e.g. - if datasetPrefix is set to 20200115_GAM_ONT_SC_rrm3
#        then the output directory will be named 20200115_GAM_ONT_SC_rrm3_8829ef1
#        where we've assumed that the commit string is 8829ef1.
#      - if the same datasetPrefix as above is used, but noCommitStr is set to 1 in the input json,
#        then the output directory will be named 20200115_GAM_ONT_SC_rrm3 i.e. without the commit string.
#        This option is useful if you update the git repo (hence changing the commit string) and want to
#        re-run parts of the pipeline on a previous run, which already has an output directory with the
#        old commit string in its name.
noCommitStr=$(jq -r '.force_no_commit_str' "$input_json")
if [ "$noCommitStr" == "1" ]; then
    datasetPrefix="$(jq -r '.dataset_prefix' "$input_json")"
else
    datasetPrefix="$(jq -r '.dataset_prefix' "$input_json")_${COMMITSTR}"
fi

# ensure that config[scratchDir] is set and set opDir
if [ "${config[scratchDir]:-""}" == "" ]; then
    echo "scratchDir is not set in config.sh" >&2
    exit 1
fi
opDir=${config[scratchDir]}/${datasetPrefix}

# set which mode the user has chosen
# -----------------------------------
useDnascentVersion4=$(jq -r '.use_dnascent_version_4 // 0' "$input_json")
useR10DoradoRefFreeModel=$(jq -r '.use_R10_dorado_ref_free_model // 0' "$input_json")
useDnascentVersion2=0

# check if both the above are set to 1,then complain and exit
if [ "$useDnascentVersion4" -eq 1 ] && [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
    echo "Error: you cannot set both use_dnascent_version_4 and use_R10_dorado_ref_free_model to 1" >&2
    exit 1
fi

# if neither of the above are set to 1, then set useDnascentVersion2 to 1
if [ "$useDnascentVersion4" -eq 0 ] && [ "$useR10DoradoRefFreeModel" -eq 0 ]; then
    useDnascentVersion2=1
fi

# fast5 and fastq directories
# ---------------------------
# NOTE: the following options are not easy to set, so please read carefully.
#       I may not have covered all cases, so try a few things if the script fails.
#       The reason I am not setting these in stone is that it is the user's responsibility to point to files suitably.
# NOTE: These options are irrelevant and must not be set if using the R10-dorado-ref-free model.
#       Scenario 1: There's one directory with all fast5 files, and that's the starting point for the pipeline,
#                   and no other steps of the pipeline have been performed i.e. only the raw currents are available.
#                   Then,
#                       * set fast5Dir to this directory
#                       * set fastQDir to "" in the input json or leave it unspecified
#                       * set fastQFile to "" in the input json or leave it unspecified
#       Scenario 2: There's a directory let's say /a/b/c with files under fast5_* folders and fastq_* folders, where
#                   "*" stands for any string i.e. raw currents available and basecalling performed. Then,
#                       * set fast5Dir to /a/b/c
#                       * set fastQDir to /a/b/c/fastq_*
#                       * set fastQFile to /a/b/c/fastq_*/*.fastq*
#                       * Turn off the basecalling step using PIPELINE_GO=0 before that step
#                       * Search for the line that sets guppySequencingSummaryFile in this script, and set that
#                         to the sequencing summary file either directly or using a suitable json input parameter,
#                         for e.g.: guppySequencingSummaryFile=/a/b/c/filename.txt
#                       * Don't bother about the guppyLogFile or the guppyConfigFile variables
#       Scenario 3: It's one of the scenarios above, but the pipeline has been partially run and stopped
#                   somewhere. Then, the options above should still work. Just set PIPELINE_GO to 0 or 1 suitably.
#                   If you've changed the code and committed it between stopping and restarting the pipeline,
#                   the commit string will have changed, which means the string opDir will have changed, which is
#                   not desirable. So change it back to the output directory you'd used earlier.
#       Scenario 4: Some other configuration of nanopore may set some other directory or file structure.
#                   Then, go through the script and set options suitably.
fast5Dir=$(jq -r '.fast5_dir // ""' "$input_json")
fastQDir=$(jq -r '.fastq_dir // ""' "$input_json")
fastQFile=$(jq -r '.fastq_file // ""' "$input_json")

# proceed further if any of the dnascent modes are chosen
if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
  # check that fast5Dir exists
  if [ ! -d "$fast5Dir" ]; then
      echo "fast5 directory $fast5Dir does not exist" >&2
      exit 1
  fi

  # if fastQDir is not set, set it to the output directory
  if [ "$fastQDir" == "" ]; then
      fastQDir=$opDir
  fi

  # if fastQFile is not set, set it to fastQDir/datasetPrefix.fastq
  if [ "$fastQFile" == "" ]; then
      fastQFile="$fastQDir"/"$datasetPrefix".fastq
  fi
else
  # check that none of these parameters are set
  if [ "$fast5Dir" != "" ] || [ "$fastQDir" != "" ] || [ "$fastQFile" != "" ]; then
      echo "fast5Dir, fastQDir, and fastQFile should not be set if not using dnascent" >&2
      exit 1
  fi
fi

# pod5 directory
# --------------
# pod5 is the input current format in the R10 dorado ref-free model.
# If you have fast5 and want to convert to pod5, refer to the run_dorado.sh script for ideas.
pod5Dir=$(jq -r '.pod5_folder // ""' "$input_json")

# check that the pod5 directory exists if the user has chosen the R10 dorado ref-free mode
if [ ! -d "$pod5Dir" ] && [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
    echo "pod5 directory $pod5Dir does not exist"
    exit 1
fi

# make a temporary directory
# --------------------------
tmpDir=$opDir/temp
mkdir -p "$tmpDir"

# names of input/output files
# ---------------------------
# If any of the DNAscent modes are used, then set some filenames
if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
  minimap2Prefix=minimap2_fastq_${datasetPrefix}
  minimap2File=$opDir/${minimap2Prefix}.sam
  minimap2BamFile=$opDir/sam_${minimap2Prefix}.bam
  minimap2SortedBamFile=$opDir/sam_${minimap2Prefix}.sorted.bam

  dnascentIndexFile=$opDir/index_${datasetPrefix}.dnascent
  dnascentDetectPrefix=$opDir/dnascent_sam_${minimap2Prefix}
  dnascentForkSensePrefix="$dnascentDetectPrefix"
  modBAMFile=${dnascentDetectPrefix}.detect.mod.sorted.bam
fi

# If the dorado mode is used, then set some filenames
if [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
  doradoPrefix=dorado_${datasetPrefix}
  dnascentDetectPrefix=$opDir/${doradoPrefix}
  dnascentForkSensePrefix="$dnascentDetectPrefix"
  modBAMFile=${dnascentDetectPrefix}.mod.sorted.bam
  modBAMFileUnfiltered=${modBAMFile/.mod.sorted.bam/.mod.unfiltered.sorted.bam}
fi

dnascentDetectFile=${dnascentDetectPrefix}.detect
modBAMForkSenseLFile=${dnascentDetectPrefix}.forkSense.mod.left.sorted.bam
modBAMForkSenseRFile=${dnascentDetectPrefix}.forkSense.mod.right.sorted.bam
dnascentForkSenseFile=${dnascentForkSensePrefix}.forkSense
dnascentForkSenseForkBdryOutputDir=$opDir/forkSenseOverallBedgraphs

guppyConfigFile=$(jq -r '.guppy_config_file // ""' "$input_json")
guppyVersion=$(jq -r '.guppy_version // "guppy_5_07"' "$input_json")
guppySequencingSummaryFile=$(jq -r '.guppy_sequencing_summary_file // ""' "$input_json")
doradoSequencingSummaryFile=$(jq -r '.dorado_sequencing_summary_file // ""' "$input_json")

# if any of the DNAscent modes are used, then set some guppy variables
# or if the dorado mode is used, then set some dorado variables
if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
  guppyLogFileDir=$opDir/guppy_log_files
  if [ "$guppyVersion" == "5.0.7" ]; then
      guppyVersion=guppy_5_07
  fi

  if [ "$guppySequencingSummaryFile" == "" ]; then
      guppySequencingSummaryFile="$opDir"/sequencing_summary.txt
  fi

  if [ "$doradoSequencingSummaryFile" != "" ]; then
      echo "dorado parameters should not be set if using the dnascent mode" >&2
      exit 1;
  fi

elif [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
  # no log file for dorado unlike guppy. I think all output is written to stdout.
  doradoSequencingSummaryFile=$(jq -r '.dorado_sequencing_summary_file // ""' "$input_json")

  if [ "$doradoSequencingSummaryFile" == "" ]; then
      doradoSequencingSummaryFile="$opDir"/sequencing_summary.txt
  fi

  if [ "$guppyConfigFile" != "" ] || [ "$guppyVersion" != "guppy_5_07" ] || [ "$guppySequencingSummaryFile" != "" ];\
  then
      echo "guppy parameters should not be set if not using the dnascent mode" >&2
      exit 1;
  fi

fi

# depending on the mode, make sure linked pipeline stages are run
# ---------------------------------------------------------------
# Some pipeline stages are linked depending on the mode i.e.
# the user must be forced to run one stage if they want to run some other stage.
# We enforce these rules in this section.

# If dnascent4 is chosen and dnascent detect is run, then the detect file
# must be converted to the dnascent 2 format for the rest of the pipeline to function.
PIPELINE_GO_DETECT=$(get_pipeline_go "$input_json" "pipeline_go_detect")
PIPELINE_GO_IF_DNASCENT4_THEN_GET_DNASCENT2_FORMAT=$(get_pipeline_go "$input_json" \
  "pipeline_go_if_dnascent4_then_get_dnascent2_format")
PIPELINE_GO_R10_REF_FREE_BASECALL_MODCALL_ALIGN=$(get_pipeline_go "$input_json" \
  "pipeline_go_R10_ref_free_basecall_modcall_align")
PIPELINE_GO_R10_REF_FREE_CONVERT_TO_DNASCENT2_FORMAT=$(get_pipeline_go "$input_json" \
  "pipeline_go_R10_ref_free_convert_to_dnascent2_format")

if [ "$useDnascentVersion4" -eq 1 ]; then

    # if dnascent detect is run, then we must get the dnascent2 format from the modBAM file
    if [ "$PIPELINE_GO_DETECT" -eq 1 ] && [ "$PIPELINE_GO_IF_DNASCENT4_THEN_GET_DNASCENT2_FORMAT" -eq 0 ]; then
        echo "In the dnascent4 workflow, you must choose convert to version 2 if you request dnascent detect v4" >&2
        exit 1
    fi

elif [ "$useR10DoradoRefFreeModel" -eq 1 ]; then

    # if the R10 dorado ref-free model is chosen, then convert to v2 must follow basecall and modcall
    if [ "$PIPELINE_GO_R10_REF_FREE_BASECALL_MODCALL_ALIGN" -eq 1 ] && \
          [ "$PIPELINE_GO_R10_REF_FREE_CONVERT_TO_DNASCENT2_FORMAT" -eq 0 ]; then
        echo "In the R10 dorado ref-free mode, modcalling data must be converted to dnascent v2" >&2 ;
        exit 1;
    fi

fi

# run the pipeline
# ================

# prepare jid variables
# ---------------------

# do a dummy job to get an initial job id
launchJob=$(sbatch --wrap="sleep 1" -p ei-short --time=00:00:10\
    -o /dev/null -e /dev/null)
jid=${launchJob##* }

# most steps in the pipeline run one after another.
# occasionally, we may want steps to run in parallel.
# so, we declare other job id variables.
# they are initialized with the dummy job id to ensure that any stage of the pipeline
# can be safely turned off without jeopardizing the whole pipeline
jid_detect_modbam="$jid"
jid_forkSense="$jid"
jid_dnascent_qc="$jid"
jid_forkSense_left_modbam="$jid"
jid_forkSense_right_modbam="$jid"
jid_nascent_read_plots="$jid"
jid_fork_sense_nascent_read_plots="$jid"

# run the pipeline
# ----------------
# set the variable below to one before the section of the pipeline
# that needs to be run. set it to zero before whichever part of the
# pipeline that needn't be run.
PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_basecall")

# basecall
if [ "$PIPELINE_GO" -eq 1 ];
then

    # check that the string r10 is present in the guppy config file parameter and that the guppy version is set to 6.5.7
    # if basecalling is requested and dnascent version 4 is chosen
    if [ "$useDnascentVersion4" -eq 1 ]; then
      if ! grep -q "r10" <<< "$guppyConfigFile"; then
          echo "a guppy config file that corresponds to an r10 model is needed if you choose a dnascent4 workflow" >&2
          exit 1
      fi
      if [ "$guppyVersion" != "6.5.7" ]; then
        echo "guppy version must be set to 6.5.7 if you choose a dnascent4 workflow" >&2
        exit 1
      else
        guppyVersion=guppy_6_57
      fi
    fi

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do plain basecalling, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    # throw an error if guppySequencingSummaryFile is not set to "$opDir"/sequencing_summary.txt
    if [ "$guppySequencingSummaryFile" != "$opDir"/sequencing_summary.txt ]; then
      {
        echo "If you want to basecall with guppy, then guppy sequencing summary file must not be set. "
        echo "We will make the file $opDir/sequencing_summary.txt for you automatically."
      } >&2
      exit 1
    fi

    exportStr=""
    exportStr="tmpDir=$tmpDir,fast5Dir=$fast5Dir,"
    exportStr+="fastQDir=$fastQDir,guppyConfigFile=$guppyConfigFile,"
    exportStr+="guppyLogFileDir=$guppyLogFileDir,fastQFile=$fastQFile,guppyVersion=${guppyVersion}"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_guppy_gpu.sh)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_align")

# perform alignment
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do plain alignment, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    exportStr=""
    exportStr="tmpDir=$tmpDir,"
    exportStr+="minimap2File=$minimap2File,"
    exportStr+="refGenomeFasta=$refGenomeFasta,"
    exportStr+="fastQFile=$fastQFile"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_minimap.sh;)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_filter_sort_before_dnascent")

# perform filtering and sorting steps before dnascent and pycoQC
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do this filtering and sorting step, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    exportStr=""
    exportStr="tmpDir=$tmpDir,"
    exportStr+="minimap2File=$minimap2File,"
    exportStr+="minimap2BamFile=$minimap2BamFile,"
    exportStr+="minimap2SortedBamFile=$minimap2SortedBamFile,"
    exportStr+="minQual=$minQual,"
    exportStr+="minLen=$minLen,"
    exportStr+="onlyPrim=$onlyPrim"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_samtools_before_dnascent.sh;)
    jid=${launchJob##* }

fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_R10_ref_free_basecall_modcall_align")

# do base and mod call if the user has chosen the R10 dorado ref-free model
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if the dorado mode is not used
    if [ "$useR10DoradoRefFreeModel" -eq 0 ]; then
        echo "To do this basecall and mod call step, you must choose the R10-dorado-ref-free mode." >&2
        exit 1
    fi

    # throw an error if doradoSequencingSummaryFile is not set to "$opDir"/sequencing_summary.txt
    if [ "$doradoSequencingSummaryFile" != "$opDir"/sequencing_summary.txt ]; then
      {
        echo "If you want to basecall with dorado, then dorado sequencing summary file must not be set. "
        echo "We will make the file $opDir/sequencing_summary.txt for you automatically."
      } >&2
      exit 1
    fi

    # add some input output files
    doradoTmpFile=$tmpDir/dorado_temp_"$(openssl rand -hex 4)".json
    jq --arg unflt "$modBAMFileUnfiltered"  --arg flt "$modBAMFile" \
      '. += {"output_unfiltered_bam_file": $unflt, "output_filtered_bam_file": $flt}' "$input_json" |\
      jq --arg seqSum "$doradoSequencingSummaryFile" \
        '.dorado_sequencing_summary_file = (.dorado_sequencing_summary_file // $seqSum)' \
        > "$doradoTmpFile"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_dorado.sh "$doradoTmpFile" ;)
    jid=${launchJob##* }
    jid_detect_modbam="$jid"
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_pycoQC")

# run pycoQC for quality control analysis
if [ "$PIPELINE_GO" -eq 1 ];
then

    # set the sequencing summary file depending on the mode used
    if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
      sequencingSummaryFile="$guppySequencingSummaryFile"
      pycoQCBamFile="$minimap2SortedBamFile"
    elif [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
      sequencingSummaryFile="$doradoSequencingSummaryFile"
      pycoQCBamFile="$modBAMFile"
    fi

    # if both a guppy and a dorado sequencing summary file are set, then we are in an undefined state, so throw an error
    if [ "${guppySequencingSummaryFile:-}" != "" ] && [ "${doradoSequencingSummaryFile:-}" != "" ]; then
      echo "Both a guppy and a dorado sequencing summary file are set. Please set only one." >&2
      exit 1
    fi

    # we cannot do a similar check for the bam file

    # set the pycoQC directory
    pycoQCDir=$opDir/pycoQCAnalysis
    mkdir -p "$pycoQCDir"

    sbatch --dependency=afterok:"$jid" \
           --mail-user="${config[email]}" \
           --mail-type=END,FAIL \
           -o "${config[job_output_logs]}"/slurm.%N.%j.out \
           -e "${config[job_output_logs]}"/slurm.%N.%j.err \
        run_pycoQC.sh \
        "$pycoQCBamFile" "$sequencingSummaryFile" \
        "$pycoQCDir" "$tmpDir";
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_R10_ref_free_convert_to_dnascent2_format")

# perform conversion to detect v2 format if the user has chosen the R10 dorado ref-free model
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if the dorado mode is not used
    if [ "$useR10DoradoRefFreeModel" -eq 0 ]; then
        echo "To do this convert to v2 step, you must choose the R10-dorado-ref-free mode." >&2
        exit 1
    fi

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_convert_modBAM_to_detect.sh "$modBAMFile" "$refGenomeFasta" "$dnascentDetectFile" ;)
    jid=${launchJob##* }

fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_make_db")

# make database
# this is an optional step. Can turn it off using PIPELINE_GO=0
if [ "$PIPELINE_GO" -eq 1 ];
then
  echo "make database step is deprecated" >&2
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_dnascent_index")

# make dnascent index
# if you are going to run dnascent, this step is mandatory
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do the dnascent index step, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    exportStr=""
    exportStr="tmpDir=$tmpDir,"
    exportStr+="fast5Dir=$fast5Dir,"
    exportStr+="dnascentIndexFile=$dnascentIndexFile,"
    exportStr+="guppySequencingSummaryFile=$guppySequencingSummaryFile,"
    exportStr+="useDnascentVersion4=$useDnascentVersion4"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_dnascent_index.sh;)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_detect")

# run dnascent detect
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do the dnascent detect step, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    # get path to hdf5 plugin needed by DNAscent
    # ------------------------------------------
    # check that you have a plugin called some_name.so here
    # if not, consult others and get the plugin
    HDF5PluginPath=$(jq -r '.hdf5_plugin_path // ""' "$input_json")

    if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
      if [ ! -f "$HDF5PluginPath/libvbz_hdf_plugin.so" ]; then
          echo "HDF5 plugin file does not exist"
          exit 1
      fi
    else
      if [ "$HDF5PluginPath" != "" ]; then
          echo "HDF5 plugin path is set but not needed"
          exit 1
      fi
    fi

    exportStr=""
    exportStr="tmpDir=$tmpDir,"
    exportStr+="minimap2SortedBamFile=$minimap2SortedBamFile,"
    exportStr+="dnascentIndexFile=$dnascentIndexFile,"
    exportStr+="refGenomeFasta=$refGenomeFasta,"
    exportStr+="dnascentDetectFile=$dnascentDetectFile,"
    exportStr+="HDF5PluginPath=$HDF5PluginPath,"
    exportStr+="minQual=$minQual,"
    exportStr+="minLen=$minLen,"
    exportStr+="useDnascentVersion4=$useDnascentVersion4"

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_dnascent_detect.sh;)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_if_dnascent4_then_get_dnascent2_format")

# convert dnascent4 format to dnascent2 format if requested
if [ "$PIPELINE_GO" -eq 1 ];
then

    # if useDnascentVersion4 is set to 0, then we print an error and exit
    if [ "$useDnascentVersion4" -eq 0 ]; then
      echo "Error: you've requested to get dnascent2 format from dnascent4, but you've set use_dnascent_version_4 to 0"\
      >&2;
      exit 1;
    fi

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        run_dnascent4_to_dnascent2.sh "$dnascentDetectFile" overwrite_and_backup;)
    jid=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_make_modbam")

# make modBAM file
if [ "$PIPELINE_GO" -eq 1 ];
then

    # throw an error if neither DNAscent modes are used
    if [ "$useDnascentVersion2" -eq 0 ] && [ "$useDnascentVersion4" -eq 0 ]; then
      echo "To do the make modBAM from detect step, you must choose one of the two dnascent modes" >&2
      exit 1
    fi

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        run_detect_to_modBAM.sh \
                        -f "$dnascentDetectFile";)
    jid_detect_modbam=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_annotate_modbam")

# set file to annotate the modBAM file, will not be run if the user has not set this stage to run
plotBedFile=$(jq -r '.plot_bed_file // ""' "$input_json")

# annotate modBAM file
if [ "$PIPELINE_GO" -eq 1 ] && [ -f "$plotBedFile" ];
then
    launchJob=$(sbatch --dependency=afterok:"$jid_detect_modbam" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       annotate_mod_bam_with_bed.sh "$modBAMFile" \
                       "${modBAMFile/mod.sorted.bam/mod.annotated.sorted.bam}" "$plotBedFile":plot_elements ;)
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_forkSense")

# run dnascent forkSense
if [ "$PIPELINE_GO" -eq 1 ];
then

    currDir=$(pwd)

    exportStr=""
    exportStr="tmpDir=$tmpDir,"
    exportStr+="dnascentDetectFile=$dnascentDetectFile,"
    exportStr+="dnascentForkSenseFile=$dnascentForkSenseFile,"
    exportStr+="scriptDir=$currDir"

    mkdir -p "$dnascentForkSenseForkBdryOutputDir";

    launchJob=$(sbatch --dependency=afterok:"$jid" \
                       -D "$dnascentForkSenseForkBdryOutputDir" \
                       --export="$exportStr" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       run_dnascent_forkSense.sh;)
    jid_forkSense=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_dnascent_qc")

# run dnascent qc
if [ "$PIPELINE_GO" -eq 1 ];
then

    # set the sequencing summary file depending on the mode used
    if [ "$useDnascentVersion2" -eq 1 ] || [ "$useDnascentVersion4" -eq 1 ]; then
      sequencingSummaryFile="$guppySequencingSummaryFile"
    elif [ "$useR10DoradoRefFreeModel" -eq 1 ]; then
      sequencingSummaryFile="$doradoSequencingSummaryFile"
    fi

    # if both a guppy and a dorado sequencing summary file are set, then we are in an undefined state, so throw an error
    if [ "${guppySequencingSummaryFile:-}" != "" ] && [ "${doradoSequencingSummaryFile:-}" != "" ]; then
      echo "Both a guppy and a dorado sequencing summary file are set. Please set only one." >&2
      exit 1
    fi

    launchJob=$(sbatch --dependency=afterok:"$jid":"$jid_detect_modbam":"$jid_forkSense" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       --mem=200G \
                       --time=23:59:59 \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                       dnascent_qc.sh "$sequencingSummaryFile" "$modBAMFile" \
                           "$dnascentForkSenseForkBdryOutputDir" "$opDir" ;)
    jid_dnascent_qc=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_make_forkSense_modbam")

# make fork sense modBAM file
if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(sbatch --dependency=afterok:"$jid_forkSense" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        run_forkSense_to_modBAM.sh \
                        -f "$dnascentForkSenseFile" \
                        -r "$refGenomeFasta" \
                        -d left ;)
    jid_forkSense_left_modbam=${launchJob##* }

    launchJob=$(sbatch --dependency=afterok:"$jid_forkSense" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        run_forkSense_to_modBAM.sh \
                        -f "$dnascentForkSenseFile" \
                        -r "$refGenomeFasta" \
                        -d right ;)
    jid_forkSense_right_modbam=${launchJob##* }
fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_plot_random_nascent_reads")

# plot a few nascent reads picked at random
if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(cd plotting_and_short_analyses;
              sbatch --dependency=afterok:"$jid_dnascent_qc" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        plot_n_nascent_reads.sh 30 "$modBAMFile" \
                            "$modBAMFile"_statistics "$opDir"/nascent_read_plots 300 0.05 ;)
    jid_nascent_read_plots=${launchJob##* }

fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_plot_random_nascent_reads_forkSense")

# plot fork sense data corresponding to the random reads from the previous step
if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(cd plotting_and_short_analyses;
              sbatch \
              --dependency=afterok:"$jid_forkSense_left_modbam":"$jid_forkSense_right_modbam":"$jid_nascent_read_plots"\
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        plot_forkSense_data_to_accompany_nascent_reads.sh "$modBAMForkSenseLFile" \
                            "$modBAMForkSenseRFile" "$opDir"/nascent_read_plots ;)
    jid_fork_sense_nascent_read_plots=${launchJob##* }

fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_plot_random_nascent_reads_collate_pdf")

# collate plots made in the previous steps along with per-read information into a pdf
if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(sbatch -p ei-short --mail-user= \
                  --dependency=afterok:"$jid_dnascent_qc":"$jid_fork_sense_nascent_read_plots"\
                  --wrap="sed '1i\detectIndex\tmean_BrdU' ${modBAMFile}_statistics > ${modBAMFile}_statistics_w_header")
    jid_temp=${launchJob##* }

    launchJob=$(cd plotting_and_short_analyses;
              sbatch \
              --dependency=afterok:"$jid_temp"\
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -p ei-medium \
                       --time=23:59:59 \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        convert_plots_to_latex.sh "$modBAMFile"_statistics_w_header \
                            "$opDir"/nascent_read_plots "$opDir"/forkSenseOverallBedgraphs \
                            "$opDir"/nascent_read_plots forkSense ;)
    jid_temp=${launchJob##* }

    launchJob=$(sbatch -p ei-short --mail-user= --dependency=afterok:"$jid_temp" \
                  --wrap="rm ${modBAMFile}_statistics_w_header")

fi

PIPELINE_GO=$(get_pipeline_go "$input_json" "pipeline_go_plot_filtered_random_nascent_reads")

# plot a few nascent reads picked at random from a filtered subset.
# note that this job depends on the dnascent qc step being finished.
# usually the filtration criterion is forks of a minimum length, alignments of a minimum length,
# and no chromosome M. Look at the command line arguments below.

if [ "$PIPELINE_GO" -eq 1 ];
then

    launchJob=$(cd plotting_and_short_analyses;
              sbatch --dependency=afterok:"$jid_dnascent_qc" \
                       --mail-user="${config[email]}" \
                       --mail-type=END,FAIL \
                       -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                       -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                        plot_filtered_n_nascent_reads_noChrM.sh 30 10 30 "$modBAMFile" \
                            "$modBAMForkSenseLFile" "$modBAMForkSenseRFile" \
                            "$dnascentForkSenseForkBdryOutputDir" \
                            "$opDir"/nascent_read_plots_forkLen10_alignLen30_noChrM ;)

fi
