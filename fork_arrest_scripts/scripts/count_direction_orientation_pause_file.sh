#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J countDirOrnPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --time=00:09:00
#SBATCH --constraint=""

# goal
# -----
# count the number of pauses in a pause file occurring on left/right forks, and on fwd/rev strands.
# optionally if BAM files are provided, we will look for nearby indels as well.

# usage
#------
# bash count_direction_orientation_pause_file.sh $pause_file [$bam_file]
# sbatch can be used in place of bash but it's not worth it as the script is very fast
# NOTE: [] denotes optional arguments, to specify them, remove the square brackets.
# pause_file: file containing pause information in our standard format, need detectIndex and (optionally) keep* columns.
#             refer to the file convert_pause_file_to_forkSense_calls.py for more information on what these columns are.
# bam_file: (optional) BAM file containing reads. if provided, we will look for nearby indels as well on the same read.

# outputs
# -------
# 14 (+2) lines: 7 numbers on 7 lines, each preceded by a line describing what the number is.
#                optionally, two more lines about indels.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -samtools -bedtools

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
pause_file=${1:-}
bam_file=${2:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 pause_file [bam_file]"
    >&2 echo "pause_file: file containing pause information in our standard format"
    >&2 echo "bam_file: (optional) BAM file containing reads"
    exit 1;
fi

# convert pause file to forkSense format. this will help us do counts.
< "$pause_file" python convert_pause_file_to_forkSense_calls.py --outputDir "$tmpDir" --discardUnkeptPauses

left_forks="$tmpDir"/leftForks_DNAscent_forkSense.bed
right_forks="$tmpDir"/rightForks_DNAscent_forkSense.bed

# count the number of pauses on left/right forks, and on fwd/rev strands
echo "Number of (kept) pauses on left forks and forward strands"
awk '{if($8 == "fwd"){print $0}}' "$left_forks" | wc -l
echo "Number of (kept) pauses on left forks and reverse strands"
awk '{if($8 == "rev"){print $0}}' "$left_forks" | wc -l
echo "Number of (kept) pauses on right forks and forward strands"
awk '{if($8 == "fwd"){print $0}}' "$right_forks" | wc -l
echo "Number of (kept) pauses on right forks and reverse strands"
awk '{if($8 == "rev"){print $0}}' "$right_forks" | wc -l
echo "Number of (kept) pauses on left forks"
< "$left_forks" wc -l
echo "Number of (kept) pauses on right forks"
< "$right_forks" wc -l
echo "Total number of (kept) pauses"
cat "$left_forks" "$right_forks" | wc -l

# if bam file exists, we will also count indels
if [ -f "$bam_file" ]; then

    echo "Total number of (kept) pauses on forks with indels"

    # form a bed file of forks
    cat "$left_forks" "$right_forks" |\
      awk -v OFS='\t' '{print $1, $2, $3, $4, "1000", ($8 == "fwd" ? "+" : "-")}' |\
        sort -k 1,1 -k2,2n > "$tmpDir"/all_forks.bed

    # extract read ids from all forks
    awk '{print $4}' "$tmpDir"/all_forks.bed | sort | uniq > "$tmpDir"/all_forks_read_ids.txt

    # make a bed file with just the first three columns of the fork bed file
    bedtools merge -i "$tmpDir"/all_forks.bed > "$tmpDir"/all_forks_3col.bed

    # get indel locations
    samtools view -N "$tmpDir"/all_forks_read_ids.txt --region-file "$tmpDir"/all_forks_3col.bed "$bam_file" |\
       python convert_bam_to_bed_but_split_by_indels.py --inThreshold 2000 --delThreshold 2000 |\
          awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' > "$tmpDir"/all_reads_split_at_indels.bed

    # form a bed file of forks and count intersections with indels
    bedtools intersect -s -a "$tmpDir"/all_forks.bed -b "$tmpDir"/all_reads_split_at_indels.bed -wa -wb |\
      awk -v OFS='\t' '$1 == $7 && $4 == $10 && $6 == $12{print $4}' |\
        sort | uniq -c | awk -v sum=0 '{if($1 > 1){sum += 1}} END{print sum}'
fi

# remove temporary directory
rm -rf "$tmpDir"