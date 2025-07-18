#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J cracDataPerGene
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Convert RNAP2 CRAC data in 10.1016/j.molcel.2022.06.021 from bigwig format to a mean-per-gene in bed file format

# usage
# -----
# sbatch convert_aiello_bigwig_to_per_gene_bed.sh $prefix_bigwig $ref_fai $gene_bed $tRNA_rDNA_bed [$chrM_optional] \
#   [$keyword_optional] [$focus_tss_or_tes_optional] [$grow_bp_optional] [$grow_which_end_optional]
# Can use bash instead of sbatch but might take a few minutes.
# [] means optional arguments. To specify an optional argument, you must specify all optional arguments to the left of
# it as well, and you must remove the square brackets for the ones you want to specify.
# prefix_bigwig: If files are at /to/file/GSM_R1_minus.bw and ...plus.bw, then prefix is /to/file/GSM_R1.
#                Aiello et al have a plus and a minus bigwig file per experimental condition.
#                So, download those, put them in a folder, and use the prefix corresponding to plus and minus of the
#                same condition here.
# ref_fai: Fasta index file for reference genome
# gene_bed: Bed file with gene annotations, must have at least six columns
# tRNA_rDNA_bed: Bed file with tRNA and rDNA annotations, must have at least six columns.
#                Used to remove data in tRNA and rDNA regions on both strands.
# chrM_optional: (default chrM) optional (regex-like) contig to exclude from analysis
#                NOTE: this is regex-like, so for e.g. chrI will match chrI, chrII, chrIII, etc.
# keyword_optional: (default "gene") optional keyword to add to the output file name. By default, the output file name
#                   is $prefix_combine_w_gene.no_tRNA_rDNA.bed. If keyword is say "ORF", then the output file name will
#                   be $prefix_combine_w_ORF.no_tRNA_rDNA.bed. This is useful when you want to analyze different
#                   annotations separately.
# focus_tss_or_tes_optional: (default "no") choices are "no", "TSS", "TES". If "TSS" or "TES" is chosen, then the
#                            gene_bed is shrunk to a zero0centered interval at the TSS or TES respectively.
#                            This is useful if you want to measure the signal only around the TSS or TES.
#                            If "no" is chosen, then the whole gene is used. NOTE that if you set this to TSS or TES,
#                            you must also set grow_bp_optional to a positive integer, otherwise all the intervals
#                            will be of size 0 and you won't capture any signal, so the program is designed to fail.
# grow_bp_optional: (default 0) grow each interval in the gene_bed file by +- so many bp. This is useful if
#                   you want to capture the signal around the element (TES or TSS) as well. As stated
#                   above, if you set focus_tss_or_tes_optional to TSS or TES, you must set this to a positive integer.
#                   If you set focus_tss_or_tes_optional to "no", then you will get an error if you set this.
# grow_which_end_optional: (default 'both') choices are 'both', '5p', '3p', 'start', 'end'.
#                           Default 'both' means each interval in the bed file is grown by +- grow_bp.
#                           If '5p' is chosen, then start -> start - grow_bp and end -> end if strand is + and
#                           start -> start and end -> end + grow_bp if strand is -.
#                           If '3p' is chosen, then start -> start and end -> end + grow_bp if strand is + and
#                           start -> start - grow_bp and end -> end if strand is -.
#                           If 'start' is chosen, then start -> start - grow_bp and end -> end irrespective of strand.
#                           If 'end' is chosen, then start -> start and end -> end + grow_bp irrespective of strand.
#                           In the above 'start' and 'end' refer to the start and end of the interval in the bed file.

# logic
# -----
# 1. Convert bigwig to bedgraph
# 2. Use bedtools unionbedg to put zero values in regions with no data.
#    Unionbedg does not work with one bedgraph so we need to give the same bedgraph twice to trick it.
# 3. Calculate mean CRAC signal per gene bed entry.
#    Refer to calculate_bed_overlap_bedgraph_signal_per_bed_entry.sh for details.
#    This is where chrM is removed as well.
# 4. Remove any gene entries that overlap with tRNAs or rDNAs on either strand.

# outputs
# -------
# One output file in bed 6+1 format with the filename ${prefix_bigwig}_combine_w_${keyword_optional}.no_tRNA_rDNA.bed
# The 6+1 format is: chr, start, end, name, score, strand, mean_per_gene
# comments start with # and no column names are included

# stop execution if any command fails
set -e

# assign arguments to variables
prefix=${1:-}
ref_fai=${2:-}
gene_bed=${3:-}
tRNA_rDNA_bed=${4:-}
chrM=${5:-chrM}
keyword=${6:-gene}
focus_tss_or_tes=${7:-no}
grow_bp=${8:-0}
grow_which_end=${9:-both}

# convert paths to absolute paths
prefix=$(realpath "$prefix")
ref_fai=$(realpath "$ref_fai")
gene_bed=$(realpath "$gene_bed")
tRNA_rDNA_bed=$(realpath "$tRNA_rDNA_bed")

# go to the main directory for scripts
cd ..;

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -python -bigWigToBedGraph -jq

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 \$prefix_bigwig \$ref_fai \$gene_bed \$tRNA_rDNA_bed \$chrM_optional \$keyword_optional "\
    "[\$focus_tss_or_tes_optional] [\$grow_bp_optional] [\$grow_which_end_optional]"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "Exiting."
    exit 1;
fi

# set up some file names
file_plus="$prefix"_plus.bw
file_minus="$prefix"_minus.bw
output_file="$prefix"_combine_w_${keyword}.no_tRNA.no_rDNA.bed

# check that input files exist
if [ ! -f "$file_plus" ]; then
    >&2 echo "ERROR: $file_plus does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$file_minus" ]; then
    >&2 echo "ERROR: $file_minus does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$ref_fai" ]; then
    >&2 echo "ERROR: $ref_fai does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$gene_bed" ]; then
    >&2 echo "ERROR: $gene_bed does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

if [ ! -f "$tRNA_rDNA_bed" ]; then
    >&2 echo "ERROR: $tRNA_rDNA_bed does not exist"
    >&2 echo "Exiting."
    exit 1;
fi

# ensure that the input bed files are valid
if [ ! "$(< "$gene_bed" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: bed file $gene_bed is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$gene_bed" python validate_bed_against_fai.py "$ref_fai" )" == "valid"  ]; then
  >&2 echo "Error: bed file $gene_bed does not have valid coordinates."
  exit 1;
fi

if [ ! "$(< "$tRNA_rDNA_bed" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: bed file $tRNA_rDNA_bed is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$tRNA_rDNA_bed" python validate_bed_against_fai.py "$ref_fai" )" == "valid"  ]; then
  >&2 echo "Error: bed file $tRNA_rDNA_bed does not have valid coordinates."
  exit 1;
fi

# check that focus_tss_or_tes is one of the allowed values
if [ "$focus_tss_or_tes" != "no" ] && [ "$focus_tss_or_tes" != "TSS" ] && [ "$focus_tss_or_tes" != "TES" ]; then
    >&2 echo "ERROR: focus_tss_or_tes must be one of 'no', 'TSS', or 'TES'"
    >&2 echo "Exiting."
    exit 1;
fi

# check that grow_bp is a number
if ! [[ "$grow_bp" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: grow_bp must be a number"
    >&2 echo "Exiting."
    exit 1;
fi

# check that if focus_tss_or_tes is set to "no", then grow_bp must be 0, and
# if focus_tss_or_tes is set to "TSS" or "TES", then grow_bp must be a positive integer
if [ "$focus_tss_or_tes" == "no" ] && [ "$grow_bp" -ne 0 ]; then
    >&2 echo "ERROR: if focus_tss_or_tes is set to 'no', then grow_bp must be 0"
    >&2 echo "Exiting."
    exit 1;
fi
if [ "$focus_tss_or_tes" != "no" ] && [ "$grow_bp" -le 0 ]; then
    >&2 echo "ERROR: if focus_tss_or_tes is set to 'TSS' or 'TES', then grow_bp must be a positive integer"
    >&2 echo "Exiting."
    exit 1;
fi

# check that grow_which_end is one of the allowed values
if [ "$grow_which_end" != "both" ] && [ "$grow_which_end" != "5p" ] && [ "$grow_which_end" != "3p" ] &&\
   [ "$grow_which_end" != "start" ] && [ "$grow_which_end" != "end" ]; then
    >&2 echo "ERROR: grow_which_end must be one of 'both', '5p', '3p', 'start', or 'end'"
    >&2 echo "Exiting."
    exit 1;
fi

# make some temp files
tmp_file_plus_bg=$(mktemp -p "$tmpDir")
tmp_file_minus_bg=$(mktemp -p "$tmpDir")
tmp_file_plus_bg_with_zero=$(mktemp -p "$tmpDir")
tmp_file_minus_bg_with_zero=$(mktemp -p "$tmpDir")

# perform the conversion from bigwig to bedgraph
bigWigToBedGraph "$file_plus" "$tmp_file_plus_bg" &
bigWigToBedGraph "$file_minus" "$tmp_file_minus_bg" &
wait

# insert zero values for any regions that have no signal.
# NOTE: we expand intervals into several intervals of 10bp each because
#       a downstream program might suffer if intervals are too large and it performs averaging that is sensitive to
#       interval size.
# NOTE: we assume the quantity measured is an intensive quantity i.e. it is the same irrespective of the interval size.
#       If the quantity were an extensive quantity, then when an interval is chopped into smaller intervals,
#       the quantity must be subdivided proportionally.
# shellcheck disable=SC2016
awk_zero_expand_cmd='
{
  if ($2 + step >= $3) {
    print $1, $2, $3, $4;
  } else {
    for (start = $2; start + step < $3; start = start + step) {
      print $1, int(start), int(start + step), $4;
    }
    print $1, int(start), int($3), $4;
  }
}'
bedtools unionbedg -i "$tmp_file_plus_bg" "$tmp_file_plus_bg" -empty -g "$ref_fai" |\
  awk -F' ' -v OFS=" " -v step=10 "$awk_zero_expand_cmd" > "$tmp_file_plus_bg_with_zero" &
bedtools unionbedg -i "$tmp_file_minus_bg" "$tmp_file_minus_bg" -empty -g "$ref_fai" |\
  awk -F' ' -v OFS=" " -v step=10 "$awk_zero_expand_cmd" > "$tmp_file_minus_bg_with_zero" &
wait

# perform bed file alterations here according to the chosen options of focus_tss_or_tes and grow_bp
tmp_gene_bed_step_1=$(mktemp -p "$tmpDir")
tmp_out_file=$(mktemp -p "$tmpDir")
# shellcheck disable=SC1010
# shellcheck disable=SC2016
mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '($6=="+" || $6=="-")' \
  then put '
    begin {
      @focus_tss_or_tes = "'"$focus_tss_or_tes"'";
    };
    if((@focus_tss_or_tes=="TSS" && $6=="+") || (@focus_tss_or_tes=="TES" && $6=="-")){
      $3=$2;
    } elif((@focus_tss_or_tes=="TSS" && $6=="-") || (@focus_tss_or_tes=="TES" && $6=="+")){
      $2=$3;
    } elif(@focus_tss_or_tes!="no" && @focus_tss_or_tes!="TSS" && @focus_tss_or_tes!="TES"){
      # we are in error if we reach here
      $3=0;
      $2=0;
    }' "$gene_bed" > "$tmp_gene_bed_step_1"

< "$tmp_gene_bed_step_1" python grow_and_split_bed_file_by_column_value_and_length.py \
  --grow-region-num-bases "$grow_bp" --fai-file "$ref_fai" --output-prefix "$tmpDir"/expand_bed \
  --allow-self-intersections --grow-which-end "$grow_which_end" > "$tmp_out_file"

# ensure that only one bed file is output
if [ "$(jq -r '. | length' "$tmp_out_file")" -ne 1 ]; then
  >&2 echo "ERROR: expected only one bed file to be output but got more than one."
  >&2 echo "Exiting."
  exit 1;
fi
tmp_gene_bed_step_2="$(jq -r '.[0].bed_file' "$tmp_out_file")"

# calculate mean signal per gene and remove any genes that overlap with tRNAs or rDNAs
bash calculate_bed_overlap_bedgraph_signal_per_bed_entry.sh "$tmp_file_plus_bg_with_zero"\
  "$tmp_file_minus_bg_with_zero" "$ref_fai" "$tmp_gene_bed_step_2" "$chrM" |\
    bedtools intersect -v -a - -b "$tRNA_rDNA_bed" |\
      insert_calling_script_header "$@" > "$output_file"

# remove temporary directory
rm -rf "$tmpDir"