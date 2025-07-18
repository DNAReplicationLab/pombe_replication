#!/usr/bin/env bash
#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J bedRegPauseReport
#SBATCH --mail-type=NONE
#SBATCH --constraint=""
#SBATCH -x t512n15
#SBATCH --time=11:59:59

# goal of program
# ================
# Print statistics associated with region(s) of interest on a reference genome,
# focusing on pause-related statistics

# program usage
# =============
# sbatch <program_name>.sh input_json_file
# can use bash in place of sbatch

# input_json_file must look like the following (fields can be in any order):
# {
#   "mod_bam": "/path/to/mod_bam_file.bam",
#   "forksense_dir": "/path/to/fork_sense_directory",
#   "pause_file": "/path/to/pause_file",
#   "bed_file": "/path/to/bed_file",
#   "op_dir": "/path/to/output_directory",
#   "dataset": "dataset1",
#   "feature": "tRNA",
#   "division": "01",
#   "analysis_label": "12jun23",
#   "mod_bam_left": "/path/to/mod_bam_left_file.bam",
#   "mod_bam_right": "/path/to/mod_bam_right_file.bam",
#   "alignment_file": "/path/to/bam_file.bam",
#   "fasta": "/path/to/fasta_file.fa", # key can also be "fasta_file" instead of "fasta"
#   "n": 10,
#   "prefix": "collection",
#   "prefix_plot_option": false,
#   "relative_direction_option": "all",
#   "delete_new_mod_bam_made_on_the_fly": false,
#   "genome_size_bp_optional": 0,
#   "AT_ratio_optional": 0,
#   "pause_sensitivity_bedgraphs_prefix_optional": "/path/to/file",
#   "restrict_fork_direction_optional": "L"
# }
#
# meanings of fields
# ==================
# NOTE: optional means you can omit the field in the json file if you don't want to use it.
# mod_bam: bam file with modification probabilities. The pipeline produces this.
# forksense_dir: directory that contains files with left forks, right forks, origins, terminations information.
#                These files are produced by forkSense in the pipeline.
# pause_file: tab-separated file with headers, containing pause information across region(s) larger than the
#             region of interest e.g.: across the whole genome.
# bed_file: coordinates on the reference genome (using half-open interval notation)
# op_dir: output directory where pdf files/plots/other files are sent.
# dataset: name of the experimental dataset
# feature: name of the feature (e.g. tRNA, rRNA, gene, etc.)
# division: label that marks if the bed_file is a subset of a larger region (e.g.: 1 for highly transcribed,
#           2 for not highly transcribed). set to "NA" if not applicable. this is used merely as a label.
# analysis_label: label that marks the analysis done on the dataset to get these pauses
# mod_bam_left: (optional) bam file with left fork probability at each thymidine, used only in per-read plotting.
# mod_bam_right: (optional) bam file with right fork probability at each thymidine, used only in per-read plotting.
# alignment_file: (optional) bam file with alignment information. This is used to plot ref vs query coordinate per read.
# fasta: (optional) fasta reference genome file. For compatibility with other scripts, you can also use
#         the key "fasta_file" instead of "fasta".
# n: (optional, default 0) number of randomly-selected reads per various categories to plot, set a small number ~10-50
# prefix: (optional, default "collection") common prefix for output file names
# prefix_plot_option: (optional, default false) if true, use common prefix for output plot folder names as well
# relative_direction_option: (optional, default "all") type of relative directions to allow,
#                             set to "all", "co-directional" or "head-on". The adjectives refer to the relative
#                             orientation between a fork and a region of interest.
# delete_new_mod_bam_made_on_the_fly: (optional, default false) if true, delete the new mod_bam file made on the fly
#                                      This is the new mod bam file formed by intersecting the input mod bam with the
#                                      input bed file. If you are making many pause reports, you will make many new
#                                      mod bam files. If you don't delete them, they will take up a lot of space.
#                                      If you are only making one pause report, you can set this to false and keep
#                                      the new mod bam file.
# genome_size_bp_optional: (optional, default 0) genome size in base pairs. Remember to include as many rDNAs as you
#                          think necessary and exclude mitochondrial genome if need be. Due to these assumptions, we
#                          do not calculate this automatically from the fasta file.
#                          By default or if set to 0, program will assume this information is not provided.
# AT_ratio_optional: (optional, default 0) AT ratio of the genome, defined as number of As and Ts divided by the
#                    total number of bases in the genome. If set to 0, the parameter is unused.
#                    If provided, the report includes the AT ratio of the region of interest (calculated from the fasta
#                    file) and the genomic AT ratio (this parameter) in the output.
#                    Downstream scripts may adjust some of our null hypotheses using
#                    this quantity as our method only calls pauses at thymidines.
#                    For more details, please refer to the script calculate_region_pause_statistics.sh.
#                    While calculating this ratio, the user must include as many rDNAs as they think necessary and
#                    must exclude the mitochondrial genome if need be, and similar such adjustments.
#                    Due to these assumptions, we do not calculate this automatically from the fasta file.
# pause_sensitivity_bedgraphs_prefix_optional: (default "" i.e. unused) prefix for pause sensitivity bedgraphs.
#                                               Let's say this prefix is set to /path/to/file, then the script
#                                               will look for the following three files: /path/to/file.all.bedgraph,
#                                               /path/to/file.left.bedgraph, /path/to/file.right.bedgraph.
#                                               These files are four-column text-separated files with no column names,
#                                               comments starting with '#' and the columns: contig, start, end, value.
#                                               The value column is the pause sensitivity value for that position.
#                                               See run_get_sgm_pause_sensitivities.sh for how to generate such files
#                                               or other details. By default, this field is not used and any associated
#                                               calculations are not performed.
# restrict_fork_direction_optional: (default "" i.e. no restriction) restrict fork direction to "L"/"R"/"lead"/"lag"
#                                    i.e. only use left, right, leading or lagging forks respectively in any
#                                    calculation. We know that in reality forks perform both leading and lagging strand
#                                    synthesis. As our data is _single-molecule_, we _can_ classify forks as leading
#                                    or lagging. Default is no restriction.
# NOTE: if mod_bam_left, mod_bam_right, fasta, and/or n are not provided, then the program will not plot reads.
#       It will still produce a pdf file with statistics and other plots, but no read plots.

# outputs
# =======
# A bunch of text files, pdf files, and optionally plots are produced in the output directory.
# As there are lots of files produced by different scripts, we are not listing them here.

# load packages
source load_package.sh -miller -latex -python -samtools -R -bedtools -jq > /dev/null 2>&1

# load configuration
source config.sh

# print node information
bash print_node_information.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# set the input json file
input_json=${1:-};

# check that at least one argument was provided
if [ -z "$input_json" ]; then
  >&2 echo "Error: no input json file provided";
  exit 1;
fi

# check that the input json file exists
if [ ! -f "$input_json" ]; then
  >&2 echo "Error: input json file does not exist";
  exit 1;
fi

# load variables from the command line
mod_bam="$(jq -r '.mod_bam' "$input_json")";
forksense_dir=$(cd "$(jq -r '.forksense_dir' "$input_json")" || exit; pwd) # set fork sense directory
pause_file="$(jq -r '.pause_file' "$input_json")";
bed_file="$(jq -r '.bed_file' "$input_json")";
op_dir_temp="$(jq -r '.op_dir' "$input_json")"
mkdir -p "$op_dir_temp"
op_dir=$(cd "$op_dir_temp" || exit; pwd) # set output directory
mod_bam_left="$(jq -r '.mod_bam_left' "$input_json")";
mod_bam_right="$(jq -r '.mod_bam_right' "$input_json")";
bam_file="$(jq -r '.alignment_file' "$input_json")";
fasta_file="$(jq -r '.fasta // .fasta_file' "$input_json")";
n_reads="$(jq -r '.n' "$input_json")";

dataset="$(jq -r '.dataset' "$input_json")";
feature="$(jq -r '.feature' "$input_json")";
division="$(jq -r '.division' "$input_json")";
analysis_label="$(jq -r '.analysis_label' "$input_json")";

# if n_reads is empty, set it to 0
if [ -z "$n_reads" ] || [ "$n_reads" == "null" ]; then
  n_reads=0;
fi

# if any of dataset, feature, division, or analysis_label are empty, set them to "NA"
if [ -z "$dataset" ] || [ "$dataset" == "null" ]; then
  dataset="NA";
fi
if [ -z "$feature" ] || [ "$feature" == "null" ]; then
  feature="NA";
fi
if [ -z "$division" ] || [ "$division" == "null" ]; then
  division="NA";
fi
if [ -z "$analysis_label" ] || [ "$analysis_label" == "null" ]; then
  analysis_label="NA";
fi

# set labels
prefix="$(jq -r '.prefix' "$input_json")";
if [ -z "$prefix" ] || [ "$prefix" == "null" ]; then
  prefix="collection";
fi

# set prefix plot option
prefix_plot_option="$(jq -r '.prefix_plot_option' "$input_json")";
if [ "$prefix_plot_option" == "true" ]; then
  prefix_plot_option="TRUE";
else
  prefix_plot_option="FALSE";
fi

# set relative direction option
relative_direction_option="$(jq -r '.relative_direction_option' "$input_json")";
if [ "$relative_direction_option" != "co-directional" ] && [ "$relative_direction_option" != "head-on" ]; then
  echo "relative_direction_option must be either all, co-directional or head-on. setting it to all" >&2;
  relative_direction_option="all";
fi

# set delete_new_mod_bam_made_on_the_fly option
delete_new_mod_bam_made_on_the_fly="$(jq -r '.delete_new_mod_bam_made_on_the_fly' "$input_json")";
if [ "$delete_new_mod_bam_made_on_the_fly" == "true" ]; then
  delete_new_mod_bam_made_on_the_fly="TRUE";
else
  delete_new_mod_bam_made_on_the_fly="FALSE";
fi

# set genome size if requested
genome_size_bp_optional="$(jq -r '.genome_size_bp_optional' "$input_json")";
if [ -z "$genome_size_bp_optional" ] || [ "$genome_size_bp_optional" == "null" ]; then
  genome_size_bp_optional=0;
fi

# set AT ratio if provided
AT_ratio_optional="$(jq -r '.AT_ratio_optional // 0' "$input_json")";

# set pause_sensitivity_bedgraphs_prefix_optional if requested
pause_sensitivity_bedgraphs_prefix_optional="$(jq -r '.pause_sensitivity_bedgraphs_prefix_optional' "$input_json")";
if [ -z "$pause_sensitivity_bedgraphs_prefix_optional" ] || \
    [ "$pause_sensitivity_bedgraphs_prefix_optional" == "null" ]; then
  pause_sensitivity_bedgraphs_prefix_optional="";
fi

# restrict based on fork direction if requested
restrict_fork_direction_optional="$(jq -r '.restrict_fork_direction_optional' "$input_json")";
if [ -z "$restrict_fork_direction_optional" ] || [ "$restrict_fork_direction_optional" == "null" ]; then
  restrict_fork_direction_optional="";
fi

# check that the input files exist
if [ ! -f "$mod_bam" ] || [ ! -d "$forksense_dir" ] || [ ! -f "$pause_file" ] || [ ! -f "$bed_file" ]; then
  >&2 echo "Error: one or more input files do not exist";
  exit 1;
fi

# check that the additional files for per-read plotting are available, or turn off the plots
read_plots_turned_off=0;
if [ ! -f "$mod_bam_left" ] || [ ! -f "$mod_bam_right" ] || [ ! -f "$fasta_file" ] || [ -z "$n_reads" ] ||\
  [ "$n_reads" -le 0 ]; then
  read_plots_turned_off=1;
fi

# if bam file is provided, check that it exists, otherwise set it to /dev/null
if [ ! -f "$bam_file" ] || [ "$bam_file" == "null" ]; then
  bam_file=/dev/null;
fi

# check that the corresponding .bai files are available for any bam files that have been provided
for bam_file in "$mod_bam" "$mod_bam_left" "$mod_bam_right" "$bam_file"; do
  if [ -f "$bam_file" ] && [ ! -f "$bam_file".bai ]; then
    >&2 echo "Error: $bam_file.bai does not exist. Please index the bam file using samtools index.";
    exit 1;
  fi
done

# check that the fasta index file exists if the plots are turned on
if [ "$read_plots_turned_off" -eq 0 ] && [ ! -f "$fasta_file".fai ]; then
  >&2 echo "Error: $fasta_file.fai does not exist. Please index the fasta file using samtools faidx.";
  exit 1;
fi

# check that the input bed file is valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --six-columns --allow-float-score --no-dot-strand)" == "valid" ];
then
  >&2 echo "Error: bed file needs at least six columns and must not have dot strand."
  exit 1;
fi

# check that the input pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# set output tex file
# the output pdf is the output tex file, but with '.tex' replaced with '.pdf'
op_tex="$op_dir"/"$prefix"_pause_report.tex

# set output mod bam file
op_mod_bam="$op_dir"/"$prefix".mod.bam

# set output pause file
op_all_pause_file="$op_dir"/"$prefix"_all_pauses # here, "all" means pauses whether we have confidence in them or not
op_no_pause_file="$op_dir"/"$prefix"_no_pauses   # forks with no pauses in them
op_elsewhere_pause_file="$op_dir"/"$prefix"_pauses_elsewhere # forks that pass through region but pause elsewhere
op_region_pause_file="$op_dir"/"$prefix"_pauses_in_region    # forks that pass through region and pause in it

# set output json report
op_region_pause_report="$op_dir"/"$prefix"_pauses_report.json    # report of statistics in this pdf in json format

# load configuration variables
source config.sh

# extract first six columns from the bed file and throw out zero-sized intervals
tmp_bed_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
< "$bed_file" grep -E -v '^browser|^track|^#' |\
   sort -k 1,1 -k2,2n |\
   awk -v OFS="\t" '{if($2!=$3){print $1,$2,$3,$4,$5,$6}}' > "$tmp_bed_file"

# if there are no regions of interest, then exit
if [ "$(head "$tmp_bed_file" | wc -l)" -eq 0 ]; then
  >&2 echo "Error: no regions of interest in the bed file";
  exit 1;
fi

# create a word count function
wc_file()
{
  cat | wc -l
}

# create functions for commonly used latex tasks
page_break()
{
    echo "\pagebreak"
}

section()
{
  echo "\section{$1}"
}

sub_section()
{
  echo "\subsection{$1}"
}

begin_verbatim()
{
  echo "\begin{verbatim}"
}

end_verbatim()
{
  echo "\end{verbatim}"
}

figure()
{
  echo "\begin{figure}[h!]"
  echo "\centering"
  echo "\includegraphics[scale=0.5]{$1}"
  echo "\end{figure}"
}

# write the latex output
{

  # make initial few pages like title, table of contents etc.
  echo "\documentclass{article}"
  echo "\usepackage{graphicx}"
  echo "\usepackage{hyperref}"
  echo "\usepackage[a4paper, margin=0.5in]{geometry}"
  # shellcheck disable=SC2028
  echo "\title{Pause report}"
  echo "\\author{${config[name]:-NA} \\\\ ${config[email]:-NA}}"
  # shellcheck disable=SC2028
  echo "\date{\today}"
  echo "\begin{document}"
  # shellcheck disable=SC2028
  echo "\maketitle"
  # shellcheck disable=SC2028
  echo "\tableofcontents"
  page_break

  section "Input json configuration file"
  begin_verbatim
  < "$input_json" fold -w 80
  end_verbatim
  page_break

  section "Source files/directories"
  begin_verbatim
  {
    echo "Mod bam file: $mod_bam";
    echo "";
    echo "Forksense directory: $forksense_dir";
    echo "";
    echo "Pause file: $pause_file";
    echo "";
    echo "Mod bam left fork probability file: $mod_bam_left";
    echo "";
    echo "Mod bam right fork probability file: $mod_bam_right";
    echo "";
    echo "Fasta file: $fasta_file";
    echo "";
    echo "Regions of interest (ROI): from $bed_file"
    echo "";
    # if bed file is not too long, then print it
    # otherwise, print a message saying that it is too long and print the first 200 lines
    if [ "$(wc -l < "$bed_file")" -le 200 ]; then
      echo "ROI bed file contents:"
      cat "$bed_file";
    else
      echo "ROI bed file contents: (too long to print, printing first 200 lines)"
      head -n 200 "$bed_file";
    fi
  } | fold -w 80
  end_verbatim
  page_break

  section "Output data files"
  echo "Relative direction option: $relative_direction_option ."
  echo "Restrict fork direction option: $restrict_fork_direction_optional"
  echo "Important note: Since a read can contain a fork in either direction, we can restrict reads in the mod bam "
  echo "file using the criteria that they overlap with the ROI and ignoring strands. But, we have no such problem "
  echo "with the pause file(s). So, the pause files below include a relative-direction between fork and region of "
  echo "interest criterion as requested by the user."
  begin_verbatim
  {
    if [ "$delete_new_mod_bam_made_on_the_fly" == "TRUE" ]; then
      echo "Mod bam file of reads overlapping with ROI (deleted after execution upon user request): N/A";
    else
      echo "Mod bam file of reads overlapping with ROI: $op_mod_bam";
    fi
    echo "";
    echo "Data from forks overlapping with ROI: $op_all_pause_file";
    echo "";
    echo "Data from forks overlapping with ROI but no pauses in them: $op_no_pause_file";
    echo "";
    echo "Data from forks overlapping with ROI but pauses not in ROI: $op_elsewhere_pause_file";
    echo "";
    echo "Data from forks overlapping with ROI and pauses in ROI: $op_region_pause_file";
    echo "";
    echo "Report of statistics in this pdf in json format: $op_region_pause_report";
    echo "";
    echo "Plots of pause statistics: $op_dir/${prefix}_pause_statistics_plots";
    echo "";
  } | fold -w 80
  end_verbatim
  page_break

  section "Count forks and pauses (at ROI and elsewhere)"
  echo "Relative direction option: $relative_direction_option"

  tmp_bed_stats=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  bash calculate_region_pause_statistics.sh -o "$tmp_bed_stats" "$op_dir" "$prefix" "$mod_bam" "$pause_file" \
             "$tmp_bed_file" "$relative_direction_option" "$genome_size_bp_optional" \
             "$pause_sensitivity_bedgraphs_prefix_optional" "$restrict_fork_direction_optional" \
             "$delete_new_mod_bam_made_on_the_fly" "$AT_ratio_optional" "$fasta_file" |\
              jq --arg d "$dataset" --arg f "$feature" --arg i "$division" --arg a "$analysis_label"\
                 --arg r "$relative_direction_option"\
                '.[] |= . + {"dataset": $d, "feature": $f, "division": $i, "analysis_label": $a, "relative_direction": $r}'|\
              tee "$op_region_pause_report" | python convert_json_to_latex.py

  # plot some pause statistics graphs
  {
    cd "plotting_and_short_analyses" || exit;
    if [ "$prefix_plot_option" == "TRUE" ]; then
      pause_stat_plot_dir="$op_dir"/"${prefix}"_pause_statistics_plots;
    else
      pause_stat_plot_dir="$op_dir"/pause_statistics_plots;
    fi
    mkdir -p "$pause_stat_plot_dir";
    Rscript plot_analysis_of_processed_pauses.R "$op_region_pause_file" "$pause_stat_plot_dir"
    cd ..;
  } > /dev/null

  # plot some meta-analyses on pauses
  {
    if [ "$prefix_plot_option" == "TRUE" ]; then
      meta_pause_stat_plot_dir="$op_dir"/"${prefix}"_pause_meta_analyses_plots;
    else
      meta_pause_stat_plot_dir="$op_dir"/pause_meta_analyses_plots;
    fi
    mkdir -p "$meta_pause_stat_plot_dir";
    bash get_collapsed_region_fork_pause_stats.sh "$op_all_pause_file" "$op_region_pause_file"\
      "$tmp_bed_file" "$meta_pause_stat_plot_dir"
  } > /dev/null

  page_break

  section "Statistics of ROI intervals"
  echo "NOTE: some of these statistics may be repeats of statistics from the section above."

  begin_verbatim

  fold -w 80 "$tmp_bed_stats"

  end_verbatim

  if [ "$read_plots_turned_off" -eq 0 ]; then

    # set output pause directories
    op_no_pause_dir="$op_dir"/no_pause
    op_elsewhere_pause_dir="$op_dir"/elsewhere_pause
    op_region_pause_dir="$op_dir"/region_pause

    if [ "$prefix_plot_option" == "TRUE" ]; then
      op_no_pause_dir="$op_dir"/"$prefix"_no_pause
      op_elsewhere_pause_dir="$op_dir"/"$prefix"_elsewhere_pause
      op_region_pause_dir="$op_dir"/"$prefix"_region_pause
    fi

    mkdir -p "$op_no_pause_dir" "$op_elsewhere_pause_dir" "$op_region_pause_dir";

    # change to plotting directory
    cd plotting_and_short_analyses || exit;

    # plot reads with no pauses
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    mlr --itsv --ocsv --ofs ' ' --skip-comments --headerless-csv-output\
        shuffle then head -n "$n_reads" then put '$read_id=splitax($detectIndex,"_")[1];'\
        then cut -o -f read_id,start,end "$op_no_pause_file"\
            > "$op_no_pause_dir"/no_pause_read_ids_subsampled

    sbatch plot_n_reads.sh "$op_no_pause_dir"/no_pause_read_ids_subsampled "$mod_bam" "$mod_bam_left"\
      "$mod_bam_right" "$forksense_dir" "$pause_file" "$fasta_file"\
      "$op_no_pause_dir" "$bam_file" > /dev/null

    # plot reads with pauses elsewhere
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    mlr --itsv --ocsv --ofs ' ' --skip-comments --headerless-csv-output\
        shuffle then head -n "$n_reads" then put '$read_id=splitax($detectIndex,"_")[1];'\
        then cut -o -f read_id,start,end "$op_elsewhere_pause_file"\
            > "$op_elsewhere_pause_dir"/elsewhere_pause_read_ids_subsampled

    sbatch plot_n_reads.sh "$op_elsewhere_pause_dir"/elsewhere_pause_read_ids_subsampled "$mod_bam" "$mod_bam_left"\
      "$mod_bam_right" "$forksense_dir" "$pause_file" "$fasta_file"\
      "$op_elsewhere_pause_dir" "$bam_file" > /dev/null

    # plot reads with pauses in the region of interest
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    mlr --itsv --ocsv --ofs ' ' --skip-comments --headerless-csv-output\
        shuffle then head -n "$n_reads" then put '$read_id=splitax($detectIndex,"_")[1];'\
        then cut -o -f read_id,start,end "$op_region_pause_file"\
            > "$op_region_pause_dir"/region_pause_read_ids_subsampled

    sbatch plot_n_reads.sh "$op_region_pause_dir"/region_pause_read_ids_subsampled "$mod_bam" "$mod_bam_left"\
      "$mod_bam_right" "$forksense_dir" "$pause_file" "$fasta_file"\
      "$op_region_pause_dir" "$bam_file" > /dev/null

    # change back to initial directory
    cd ..;

    section "List of folders with other related files/documents"
    begin_verbatim
    {
      echo "Reads with at least one fork with no pause"
      echo "$op_no_pause_dir"/no_pause_read_ids_subsampled
      echo ""
      echo "Reads with at least one fork with pause but not in ROI"
      echo "$op_elsewhere_pause_dir"/elsewhere_pause_read_ids_subsampled
      echo ""
      echo "Reads with at least one fork with pause and at ROI"
      echo "$op_region_pause_dir"/region_pause_read_ids_subsampled
    }  | fold -w 80
    end_verbatim

  fi

  echo "\end{document}"

} > "$op_tex"

# convert to pdf
# pdflatex sometimes needs to be run twice, so we are going to run it twice every time. This is not a bug.
pdflatex -output-directory="$op_dir" "$op_tex"
pdflatex -output-directory="$op_dir" "$op_tex"

# remove temporary files
rm "$tmp_bed_file" "$tmp_bed_stats"