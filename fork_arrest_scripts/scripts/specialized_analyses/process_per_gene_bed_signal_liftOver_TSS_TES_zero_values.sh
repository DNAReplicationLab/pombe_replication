#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J processPerGeneBedVarious
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# Goal
# -----
# Process per gene bed signal in various ways for input into our pause correlation pipeline
# (this is the pipeline where pauses already found and other datasets are correlated with them)

# Logic
# ------
# - Do the following two times, once for the gene annotation and once for the transcript annotation (which we call ORF-Ts)
# -- first, convert a mean signal into a sum signal by multiplying by the gene length
# -- extract just the TSS and TES
# -- liftOver the signal to a new genome
# - remove zero values from the ORF-Ts file.
# - gather zero values from the ORF-Ts file and the gene annotation file into a new file and grow intervals.

# usage
#------
# bash process_per_gene_bed_signal_liftOver_TSS_TES_zero_values.sh $gene_bed $orf_ts_bed $chain_file \
#   $fasta_fai $output_dir [$interval_size] [$bed_type]
# NOTE: [] indicates optional arguments
# - gene_bed: bed file containing a mean signal per gene, probably created using our
#             convert_aiello_bigwig_to_per_gene_bed.sh script. It is important that the signal is a mean signal
#             i.e. normalized by the length of the gene.
# - orf_ts_bed: same as above but mean signal per transcript. (Gene and transcripts are different annotations)
# - chain_file: chain file for liftOver.
# - fasta_fai: fasta fai file for the new genome
# - output_dir: directory to save the output files
# - interval_size: (default 250bp) size of the interval to grow the zero values in the gene annotation file
# - bed_type: (default "whole") the two bed files we use (gene_bed, orf_ts_bed) could be the whole gene or just areas
#             centered at the TSS or the TES. This parameter specifies which one it is. It can be "whole",
#             "TSS_centered", or "TES_centered".

# outputs
# -------

# stop execution if any command fails
set -e

# get current working directory and go to parent directory
cwd=$(pwd)
cd ..

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -miller -liftOver

# check that the liftOver binary is available
if [ -z "${liftOver:-}" ]; then
  >&2 echo "ERROR: liftOver binary not found"
  exit 1;
fi

# load configuration
source config.sh

# change back to the current working directory
cd "$cwd"

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
gene_bed=$1
orf_ts_bed=$2
chain_file=$3
fasta_fai=$4
output_dir=$5
interval_size=${6:-250}
bed_type=${7:-"whole"}

# check that the correct number of arguments were provided
if [ "$#" -lt 5 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash process_per_gene_bed_signal_liftOver_TSS_TES_zero_values.sh \$gene_bed \$orf_ts_bed "\
    "\$chain_file \$fasta_fai \$output_dir [\$interval_size] [\$bed_type]"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "[\$interval_size] is optional and defaults to 250"
    >&2 echo "[\$bed_type] is optional and defaults to 'whole'"
    exit 1;
fi

# check that the first four arguments are files that exist
for file in "$gene_bed" "$orf_ts_bed" "$chain_file" "$fasta_fai"; do
    if [ ! -f "$file" ]; then
        >&2 echo "ERROR: file $file does not exist"
        exit 1;
    fi
done

# check that bed_type is one of the allowed values
if [ "$bed_type" != "whole" ] && [ "$bed_type" != "TSS_centered" ] && [ "$bed_type" != "TES_centered" ]; then
  >&2 echo "ERROR: bed_type must be one of 'whole', 'TSS_centered', or 'TES_centered'"
  exit 1;
fi

# make output directory if it doesn't exist
mkdir -p "$output_dir"

# make a list of new files that are made in this script
list_of_new_files=$(mktemp --tmpdir="$tmpDir")

for bed_file in "$gene_bed" "$orf_ts_bed"; do
  for suffix in "TSS" "TES"; do

    # skip the appropriate suffix if the bed type is not "whole"
    if [ "$bed_type" != "whole" ] && [ "$suffix"_centered != "$bed_type" ]; then
      continue
    fi

    tmp_file_0=$(mktemp --tmpdir="$tmpDir")
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    mlr --tsv --implicit-csv-header --skip-comments --headerless-csv-output filter '($6=="+" || $6=="-")' \
      then put '
        begin {
          @suffix = "'"$suffix"'";
          @bed_type = "'"$bed_type"'";
        };
        $7 = $7*($3-$2);
        if(@bed_type=="whole"){
          if(($6=="+" && @suffix=="TSS") || ($6=="-" && @suffix=="TES")){
            $3=$2;
          } elif(($6=="+" && @suffix=="TES") || ($6=="-" && @suffix=="TSS")){
            $2=$3;
          }
        } elif((@suffix=="TSS" && @bed_type=="TSS_centered") || (@suffix=="TES" && @bed_type=="TES_centered")){
          bed_center = int(($2+$3)/2);
          $2=bed_center;
          $3=bed_center;
        }' "$bed_file" > "$tmp_file_0"

    new_bed="$output_dir"/$(basename "$bed_file" .bed).sumCRAC."$suffix".zero_size.bed
    unmapped=${new_bed/zero_size.bed/zero_size.unmapped.bed}

    # liftOver the bed file to the new genome
    $liftOver "$tmp_file_0" "$chain_file" "$new_bed" "$unmapped"

    # make a new file with only non-zero values
    new_bed_no_zero=${new_bed/zero_size.bed/zero_size.no_zero_values.bed}
    awk '$7!=0' "$new_bed" > "$new_bed_no_zero"

    # add the new file to the list of new files
    {
      echo "$new_bed";
      echo "$new_bed_no_zero";
      echo "$unmapped";
    } >> "$list_of_new_files"

  done

done

# make a zero value bed file
for suffix in "TSS" "TES"; do

  # skip the appropriate suffix if the bed type is not "whole"
  if [ "$bed_type" != "whole" ] && [ "$suffix"_centered != "$bed_type" ]; then
    continue
  fi

  tmp_file_1=$(mktemp --tmpdir="$tmpDir")
  tmp_file_2=$(mktemp --tmpdir="$tmpDir")

  new_gene_bed="$output_dir"/$(basename "$gene_bed" .bed).sumCRAC."$suffix".zero_size.bed
  new_orf_ts_bed="$output_dir"/$(basename "$orf_ts_bed" .bed).sumCRAC."$suffix".zero_size.bed
  new_orf_ts_bed_no_zero=${new_orf_ts_bed/zero_size.bed/zero_size.no_zero_values.bed}

  {
    awk '$7==0' "$new_gene_bed"
    awk '$7==0' "$new_orf_ts_bed"
  } |  bedtools slop -i - -b "$interval_size" -g "$fasta_fai"  > "$tmp_file_1"

  bedtools slop -i "$new_orf_ts_bed_no_zero" -b "$interval_size" -g "$fasta_fai" > "$tmp_file_2"
  output_bed=${new_gene_bed/zero_size.bed/zero_values.bed}
  output_bed=${output_bed/gene/gene_OR_ORF}

  bedtools subtract -s -a "$tmp_file_1" -b "$tmp_file_2" > "$output_bed"

  echo "$output_bed" >> "$list_of_new_files"
done

# for each new file, insert the calling script header
while IFS= read -r file; do
  < "$file" insert_calling_script_header "$@" > "$file".tmp
  mv "$file".tmp "$file"
done < "$list_of_new_files"

# remove temporary directory
rm -rf "$tmpDir"
