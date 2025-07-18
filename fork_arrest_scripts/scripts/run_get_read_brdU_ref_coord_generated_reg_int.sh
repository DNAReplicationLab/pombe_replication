#!/bin/bash

#SBATCH --mem-per-cpu=7G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J aggModBamRegInt
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1:59:59

# goal
# ====
# Generate a list of non-overlapping genomic windows of a specified size across entire genome, and calculate
# mean analogue density per read id for each window if any part of the read overlaps with the window.

# usage
# =====
# sbatch run_get_read_brdU_ref_coord_generated_reg_int.sh mod_bam fasta_file fasta_index_file window_size\
#    output_file temp_directory
# NOTE: \ means that the command continues on the next line
# NOTE: might be able to use bash in place of sbatch.
# mod_bam: modBAM file with reads mapped to reference containing analogue probability per base of interest per read
# fasta_file: fasta file with reference genome e.g. sacCer3.fa
# fasta_index_file: fasta index file with reference genome e.g.: sacCer3.fa.fai
# window_size: size of regular genomic windows to be created (in base pairs).
#   NOTE: windows at the end of a contig may be smaller.
# output_file: output file with mean analogue density per read id for each window
#         Output format is tab-separated with no headers.
#         Columns are contig, start, end, A_C_G, T, read_id, analogue_mean, num_bases.
#         The first three columns are window coordinates.
#         Fourth is a string with the number of A, C, and G bases in the window separated by underscores.
#         Fifth is the number of T bases in the window.
#         Sixth is the read id.
#         Seventh is the mean analogue density over the bases of the read in the window
#         Eighth is the total number of mod and unmod bases in the genomic window.
# temp_directory: temporary directory to store intermediate files

# stop execution if any step fails
set -e

# get config directory info
source config.sh

# read command line inputs and set variables
if [ "$#" -ne 6 ]; then
    echo "Usage: sbatch <program_name.sh> mod_bam fasta_file fasta_index_file window_size output_file temp_directory"
    exit;
fi

# set input and output directories, files
modBAM=$1
fastaFile=$2
fastaIndexFile=$3
windowSize=$4
outputFile=$5

# check that input files exist
if [ ! -f "$modBAM" ]; then
    echo "Error: modBAM file $modBAM does not exist."
    exit;
fi

if [ ! -f "$modBAM".bai ]; then
    echo "Error: modBAM file index $modBAM.bai does not exist."
    exit;
fi

if [ ! -f "$fastaFile" ]; then
    echo "Error: fasta file $fastaFile does not exist."
    exit;
fi

if [ ! -f "$fastaIndexFile" ]; then
    echo "Error: fasta index file $fastaIndexFile does not exist."
    exit;
fi

# check that window size is a positive integer
if ! [[ "$windowSize" =~ ^[0-9]+$ ]]; then
    echo "Error: window size $windowSize is not a positive integer."
    exit;
fi

# make temporary directory if need be
tmpDir=${6:-"${config[scratchDir]:-}"/tmp}
mkdir -p "$tmpDir";
tmpDir=$(realpath "$tmpDir")

# make two temporary bed files
tmpBedFile=$(mktemp "$tmpDir"/tmp_bed_file_XXXXXX.bed)
tmpBedFileComp=$(mktemp "$tmpDir"/tmp_bed_file_comp_XXXXXX.bed)

# load bedtools and seqtk
source load_package.sh -bedtools -seqtk

# make a bed file with windows of specified size from the fasta index file
bedtools makewindows -g "$fastaIndexFile" -w "$windowSize" > "$tmpBedFile"

# use seqtk to get the composition of each window
seqtk comp -r "$tmpBedFile" "$fastaFile" | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$5"_"$6"\t"$7}' > "$tmpBedFileComp"

# run run_get_read_brdU_ref_coord_int.sh, the job that performs averaging per read per window
launchJob=$(sbatch --mail-user="${config[email]}" \
                   --mail-type=END,FAIL \
                   -o "${config[job_output_logs]}"/slurm.%N.%j.out \
                   -e "${config[job_output_logs]}"/slurm.%N.%j.err \
                   run_get_read_brdU_ref_coord_int.sh "$modBAM" "$tmpBedFileComp" "$outputFile" "$tmpDir" ;)

jid=${launchJob##* }

# delete temporary files after job is done
sbatch --dependency=afterok:"$jid" -p ei-short --wrap="rm $tmpBedFile $tmpBedFileComp"