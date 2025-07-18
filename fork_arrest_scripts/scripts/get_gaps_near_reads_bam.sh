#!/bin/bash

#SBATCH --mem-per-cpu=100G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J getGapsNearReadsBam
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# goal
# -----
# Mark where indels exist in BAM files (by indels we mean insertions or deletions
# above a certain length in single molecules).

# usage
#------
# bash get_gaps_near_reads_bam.sh $bam_file [$size_threshold] [$region]
# NOTE: [] denotes optional arguments. To specify them, remove the brackets and you must also
# specify the arguments that come before them.
# bam_file: BAM file containing reads and alignment information.
#           Important: if you are passing mod BAM files, you have to pass files with alignment information,
#           and not just the mod calls in ref-anchored files because the CIGAR strings in these may not
#           contain alignment information.
# size_threshold: indels called only above this threshold. Default is 2000 bp.
# region: restrict analysis to reads passing through this region, in format "contig" or "contig:start-end".
#         Default is the whole genome. NOTE: reads with any overlap with this region are considered.

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
source load_package.sh -python -samtools

# load configuration
source config.sh

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
bam_file=${1:-}
size_threshold=${2:-2000}
region=${3:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash get_gaps_near_pauses.sh bam_file [size_threshold] [region]"
    >&2 echo "[] means optional arguments. To use, remove brackets and specify all arguments that come before them."
    >&2 echo "bam_file: BAM file containing reads and alignments"
    >&2 echo "size_threshold: indels called only above this threshold. Default is 2000 bp."
    >&2 echo "region: restrict to region, in format 'contig' or 'contig:start-end'. Default is the whole genome."
    exit 1;
fi

# check that the input file exists
if [ ! -f "$bam_file" ]; then
    >&2 echo "ERROR: BAM file not found: $bam_file"
    exit 1;
fi

# check that region is of the correct format
if [ -n "$region" ]; then
    if [[ ! "$region" =~ ^[A-Za-z0-9_.]+:[0-9]+-[0-9]+$ ]] && [[ ! "$region" =~ ^[A-Za-z0-9_.]+$ ]]; then
        >&2 echo "ERROR: region is not of the correct format: $region"
        >&2 echo "Correct format is 'contig:start-end'"
        >&2 echo "NOTE: if your contig contains unusual characters, the script may not work"
        exit 1;
    fi
fi

# get indel locations
# and connect the 5' and 3' ends of only the intermediary pieces
# shellcheck disable=SC2086
samtools view "$bam_file" $region |\
   python convert_bam_to_bed_but_split_by_indels.py --inThreshold "$size_threshold" --delThreshold "$size_threshold" |\
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
            }'  | insert_calling_script_header "$@"