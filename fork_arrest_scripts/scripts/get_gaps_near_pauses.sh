#!/bin/bash

#SBATCH --mem-per-cpu=100G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J getGapsNearPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --time=11:59:59

# goal
# -----
# Given a pause file, extract forks, and mark where indels exist in these forks (by indels we mean insertions or
# deletions above a certain length in single molecules).

# usage
#------
# bash get_gaps_near_pauses.sh $pause_file $bam_file
# pause_file: file containing pause information in our standard format, need detectIndex and (optionally) keep* columns.
#             refer to the file convert_pause_file_to_forkSense_calls.py for more information on what these columns are.
# bam_file: BAM file containing reads and alignment information.

# outputs
# -------
# bed file with columns of contig, start, end, read id to standard output.
# end - start = 1 bp for insertions and much longer than 1 bp for deletions.
# e.g.
# chrIII  72684   72685   c83b1826-db6a-4d3c-9808-5f41895ae6c0
# chrIII  72685   72686   0556bfe1-39ac-4b9b-9b61-ebb134839657
# chrIII  72686   72687   1fe81127-730e-4674-b8a6-4e87ae2d01e7
# chrIII  72665   72666   07ec3f69-f66f-42ac-90d6-5ef70aa227f9
# chrXIV  423092  436014  b57accf9-9bbf-4e05-9b1d-49a49221ed50
# chrXIV  425380  438313  15709b43-a5b9-4469-a2fa-4df7b3b939a8
# NOTE: the above is not real data
# In the above, there are insertions around 72680, whereas the last two rows are deletions.

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
if [ "$#" -lt 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash get_gaps_near_pauses.sh pause_file bam_file"
    >&2 echo "pause_file: file containing pause information in our standard format"
    >&2 echo "bam_file: BAM file containing reads and alignments"
    exit 1;
fi

# check that the two input files exist
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file not found: $pause_file"
    exit 1;
fi
if [ ! -f "$bam_file" ]; then
    >&2 echo "ERROR: BAM file not found: $bam_file"
    exit 1;
fi

# convert pause file to forkSense format. this will help us do counts.
< "$pause_file" python convert_pause_file_to_forkSense_calls.py --outputDir "$tmpDir" --discardUnkeptPauses

left_forks="$tmpDir"/leftForks_DNAscent_forkSense.bed
right_forks="$tmpDir"/rightForks_DNAscent_forkSense.bed

# form a bed file of forks
cat "$left_forks" "$right_forks" |\
  awk -v OFS='\t' '{print $1, $2, $3, $4}' |\
    sort -k 1,1 -k2,2n > "$tmpDir"/all_forks.bed

# extract read ids from all forks
awk '{print $4}' "$tmpDir"/all_forks.bed | sort | uniq > "$tmpDir"/all_forks_read_ids.txt

# make a bed file with just the first three columns of the fork bed file
bedtools merge -i "$tmpDir"/all_forks.bed > "$tmpDir"/all_forks_3col.bed

# get indel locations
# and connect the 5' and 3' ends of only the intermediary pieces
samtools view -N "$tmpDir"/all_forks_read_ids.txt --region-file "$tmpDir"/all_forks_3col.bed "$bam_file" |\
   python convert_bam_to_bed_but_split_by_indels.py --inThreshold 2000 --delThreshold 2000 |\
      grep -E -v '^#|^browser|^track' |\
        awk -v OFS='\t' '{print $2, $3, $1, $4}' |\
          uniq -f 2 --all-repeated=separate  |\
            awk -v OFS='\t' '
            {
              if ($1 != "" && prev_start != "" && $2 != "" && prev_end != "") {
                if(prev_end < $1){
                  print $3, prev_end, $1, $4
                } else {
                  print $3, prev_end, $1 + 1, $4
                }
              }
              prev_start = $1
              prev_end = $2
            }' | insert_calling_script_header "$@"

# remove temporary directory
rm -rf "$tmpDir"