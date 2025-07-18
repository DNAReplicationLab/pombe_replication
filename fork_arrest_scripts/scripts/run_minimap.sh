#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 16
#SBATCH -p ei-medium
#SBATCH -J mnmp2Full
#SBATCH --mail-type=END,FAIL
#SBATCH --time=09:59:59

set -e

# load minimap2 and git details
source load_package.sh -minimap2pt24
source load_git_repo_labels.sh

nThreads=16

{
  # add commit string
  echo -e "@CO\tfrom commit ${COMMITSTR}";
  minimap2 -L -a -x map-ont -t $nThreads ${refGenomeFasta} ${fastQFile};
} > ${minimap2File}