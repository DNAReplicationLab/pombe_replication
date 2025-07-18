#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J markHighPauseSitesRevLift
#SBATCH --mail-type=END,FAIL
#SBATCH --time=11:59:59

# Goal
# -----
# Find sites with high pause number relative to expectation, and lift them over from W303 to SacCer3.

# usage
#------
# bash mark_high_pause_sites_and_reverse_liftOver.sh $pause_file $sensitivity_file $chain_file $fasta_fai \
#    $output_directory [$sig]
# NOTE: [] means optional argument
# NOTE: \ is used to indicate that the command continues on the next line.
# - pause_file: pause file in the usual TSV format produced by our pause detection pipeline.
#               See convert_pause_file_to_bed.py for further details on the format.
# - sensitivity_file: sensitivity file in the usual bedgraph format produced by our pipeline.
#                     We want the "all" sensitivity file i.e. the null hypothesis for all forks.
#                     This file is usually a signal over regular windows of size 10 bp or so.
# - chain_file: chain file for liftOver.
# - fasta_fai: fasta index file for the genome.
# - output_dir: directory to save the output files
# - sig: (optional) significance threshold for high pause sites. Default is 3.0 i.e. the high pause sites are those with
#        pause number 3 standard deviations above the expectation.

# outputs
# -------
# - One bed file with ratio of (observed-expected)/sqrt(expected) in windows equal to the sensitivity file.
#   (if expected = 0 somewhere, we will not output that window). This bedfile will have overlapping windows.
# - One bed file (in BED3 format) with list of sites with high pause number relative to expectation in SacCer3.
# - One bed file (in BED3 format) with list of sites with high pause number relative to expectation in W303.
# - One unmapped file with sites that could not be lifted over.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -bedtools -liftOver -python

# check that the liftOver binary is available
if [ -z "${liftOver:-}" ]; then
  >&2 echo "ERROR: liftOver binary not found"
  exit 1;
fi

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

# assign arguments to variables
pause_file=$1
sensitivity_file=$2
chain_file=$3
fasta_fai=$4
output_dir=$5
sig=${6:-3}

# check that the correct number of arguments were provided
if [ "$#" -lt 5 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash mark_high_pause_sites_and_reverse_liftOver.sh \$pause_file \$sensitivity_file \$chain_file "\
"\$fasta_fai \$output_directory [\$sig]"
    >&2 echo "NOTE: [] means optional argument"
    >&2 echo "pause_file: pause file in the usual TSV format produced by our pause detection pipeline."
    >&2 echo "sensitivity_file: sensitivity file in the usual bedgraph format produced by our pipeline."
    >&2 echo "chain_file: chain file for liftOver."
    >&2 echo "fasta_fai: fasta index file for the genome."
    >&2 echo "output_dir: directory to save the output files"
    >&2 echo "sig: (optional) significance threshold for high pause sites. Default is 3.0."
    >&2 echo "For more details, see the script header."
    exit 1;
fi

# check that the first three arguments are files that exist
for file in "$pause_file" "$sensitivity_file" "$chain_file" "$fasta_fai"; do
    if [ ! -f "$file" ]; then
        >&2 echo "ERROR: file $file does not exist"
        exit 1;
    fi
done

# make output directory if it doesn't exist
mkdir -p "$output_dir"
output_dir=$(realpath "$output_dir")

# check that the input pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# convert pause file to bed 3 format
pause_bed="$tmpDir"/pause.bed
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
  grep -E -v "^browser|^track|^#" | awk -v OFS="\t" '{print $1, $2, $3, $4}' |\
    bedtools sort -i - -faidx "$fasta_fai" > "$pause_bed"

# convert sensitivity bedgraph to bed 3+1 format
sensitivity_bed="$tmpDir"/sensitivity.bed
< "$sensitivity_file" grep -v -E '^browser|^track|^#' |\
  awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4}' |\
    bedtools sort -i - -faidx "$fasta_fai" > "$sensitivity_bed"

# check that the sensitivity file is valid
validity=$(< "$sensitivity_bed" awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' |\
  python validate_bed_against_fai.py --check-genome-cov-to-tol 100  --check-no-overlap "$fasta_fai")
if [ "$validity" != "valid" ]; then
  >&2 echo "ERROR: sensitivity file $sensitivity_file is not valid"
  >&2 echo "This may not mean that the file is an invalid bedgraph, but it does not meet the requirements "
  >&2 echo "for this script."
  exit 1;
fi

if [ ! "$(< "$sensitivity_bed" awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' | python validate_bed_format.py )" == "valid" ]; then
  >&2 echo "Error: sensitivity file appears to be invalid!"
  exit 1;
fi

n_lines_not_10_bp=$(< "$sensitivity_bed" awk '$3-$2!=10' | wc -l)
n_lines_fasta_fai=$(< "$fasta_fai" wc -l)

if [ "$n_lines_not_10_bp" -gt  "$n_lines_fasta_fai" ]; then
  >&2 echo "Error: sensitivity has mis-sized windows!"
  exit 1;
fi

echo "Finished checks. Proceeding with analysis."

# - use bedtools make windows to generate a bed file of windows
# - sum the sensitivity file over the windows (cg win means coarse-grained window)
# - count number of pauses over windows
# - calculate enrichment i.e. (observed - expected)/sqrt(expected)
pauses_sensitivity_cg_win="$output_dir"/pauses_sensitivity_cg_win.bed
bedtools makewindows -g "$fasta_fai" -w 1000 -s 10 |\
  bedtools map -a - -b "$sensitivity_bed" -c 4 -o sum |\
    bedtools intersect -a - -b "$pause_bed" -c |\
      awk -v OFS="\t" '{if($4 > 0){print $1, $2, $3, $4, $5, ($5 - $4)/sqrt($4)} else {print $1, $2, $3, $4, $5, "NA"}}'|\
        insert_calling_script_header "$@" > "$pauses_sensitivity_cg_win"

# now, gather sites with high pause number relative to expectation.
# also include read ids for these sites
high_pause_sites="$output_dir"/high_pause_sites.bed
grep -v -E '^browser|^track|^#' "$pauses_sensitivity_cg_win" |\
  awk -v sig="$sig" -v OFS="\t" '$6 >= sig{print $1, $2, $3}' |\
    bedtools sort -i - -faidx "$fasta_fai" |\
      bedtools merge -d 100 -i - |\
        bedtools map -a - -b "$pause_bed" -c 4 -o distinct |\
          grep ',' |\
            sed '1i#NOTE: we are rejecting regions with fewer than two pauses' |\
              insert_calling_script_header "$@" > "$high_pause_sites"

# set names for the high pause sites
high_pause_sites_center="$output_dir"/high_pause_sites.named.bed
awk -v OFS="\t" -v a=0 '!/^#/{a+=1;print $1, $2, $3, "site_" a "_" $1 "_" $2 "_" $3 "_w303"}' \
  "$high_pause_sites" > "$high_pause_sites_center"

# now liftover
new_high_pause_sites_center="$output_dir"/high_pause_sites.named.sacCer3.bed
unmapped_center="$output_dir"/high_pause_sites.named.sacCer3.unmapped.bed

# liftOver the bed file to the new genome
$liftOver -minMatch=0.5 "$high_pause_sites_center" "$chain_file" "$new_high_pause_sites_center" "$unmapped_center"

# sort the new high pause sites
bedtools sort -i "$new_high_pause_sites_center" -faidx "$fasta_fai"|\
  insert_calling_script_header "$@" > "$new_high_pause_sites_center".sorted
mv "$new_high_pause_sites_center".sorted "$new_high_pause_sites_center"

# remove temporary directory
rm -rf "$tmpDir"