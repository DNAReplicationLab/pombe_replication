#!/bin/bash
#SBATCH --mem=50G
#SBATCH -c 16
#SBATCH -p ei-medium
#SBATCH -J sambfnas
#SBATCH --mail-type=END,FAIL
#SBATCH --time=4:00:00

set -e

# load samtools
source load_package.sh -samtools

# samtools filtering step - subsample, and ensure
# only reads of a min alignment quality and
# whose reference sequence is above a min length
# are let through.

if [ $onlyPrim -eq 1 ];
then
  # get flags ready for primary read filtering so only primary reads are allowed through
  secFlag=256
  supplFlag=2048
  invPrimFlag=$((secFlag + supplFlag))
else
  # no flag-based filtering
  invPrimFlag=00
fi

samtools view -e "rlen>=$minLen&&mapq>=$minQual" -Sb -F $invPrimFlag -o $minimap2BamFile $minimap2File;

samtools sort -@ 16 -o $minimap2SortedBamFile -T $tmpDir $minimap2BamFile
# -T temp means some temp files are created with the specified prefix of 'temp'

samtools index $minimap2SortedBamFile
# above step produces sam_blah.bam.bai