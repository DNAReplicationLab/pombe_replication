#!/bin/bash
#SBATCH --mem=20G
#SBATCH -p ei-medium
#SBATCH -J convertModBamToDetect
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59

# goal
# ----
# convert modBAM file to detect format

# input/output
# ------------
# 1. modBAM file
# 2. linear reference genome in fasta format
# 3. output file name, output will be in .detect format

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch run_convert_modBAM_to_detect.sh \$modBAMFile \$referenceGenome \$outputFile" >&2
    echo -e "\t* Can also use sbatch in place of bash."
    echo -e "\t* modBAMFile: modBAM file to convert to detect format."
    echo -e "\t* referenceGenome: linear reference genome in fasta format."
    echo -e "\t* outputFile: output file name, output will be in .detect format."
    exit 1;
fi

modBAMFile="$1"
referenceGenome="$2"
outputFile="$3"

# check that the input files exist
if [ ! -f "$modBAMFile" ]; then
    echo "modBAM file $modBAMFile does not exist." >&2
    exit 1;
fi

if [ ! -f "$referenceGenome" ] || [ ! -f "$referenceGenome.fai" ]; then
    echo "Reference genome $referenceGenome (and/or the .fai file) does not exist." >&2
    exit 1;
fi

# make any output directories if they don't exist
mkdir -p "$(dirname "$outputFile")"

# load samtools and python
source load_package.sh -python -samtools

# convert modBAM to detect
samtools view -h "$modBAMFile" |\
  python convert_modBAM_to_detect.py --fasta "$referenceGenome" |\
    python block_mismatched_detect_lines.py > "$outputFile"