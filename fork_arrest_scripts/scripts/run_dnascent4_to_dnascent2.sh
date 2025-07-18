#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J runDnascent4to2
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# convert a detect file from dnascent v4 format to v2 using only the BrdU signal

# usage
#------
# bash run_dnascent4_to_dnascent2.sh $inputFile [mode]
# (can use sbatch or bash to run)
# where:
# $inputFile: input detect file in dnascent v4 format
# mode: (default overwrite_and_backup) in overwrite_and_backup mode, the v4 detect file is replaced by the v2 detect
#        file and the original v4 detect file is backed up with a .bkup extension. As of now, only the default
#        mode is supported. So, if you specify any other mode, the script will throw an error.
# NOTE that [] denotes optional arguments. If you want to specify such an argument, you have to get rid
#      of the [], and specify the argument.

# outputs
# -------
# a detect file in dnascent v2 format. In the default mode, the v4 detect file is replaced by the v2 detect file,
# and the original v4 detect file is backed up with a .bkup extension.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

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
input_file=${1:-}
mode=${2:-overwrite_and_backup}

# set an output file
random_str=$(openssl rand -hex 8)
output_file=$input_file."$random_str"

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash run_dnascent4_to_dnascent2.sh \$inputFile [mode]"
    >&2 echo "where:"
    >&2 echo "\$inputFile: input file in dnascent v4 format"
    >&2 echo "mode: (default overwrite_and_backup) in default mode, the v4 detect file is replaced by the "
    >&2 echo "       v2 detect file and the original v4 detect file is backed up with a .bkup extension."
    >&2 echo "       As of now, only the default mode is supported."
    >&2 echo "can use sbatch or bash to run"
    >&2 echo "NOTE that [] denotes optional arguments. If you want to specify such an argument, you have to get rid"
    >&2 echo "      of the [], and specify the argument."
    exit 1;
fi

# check that the input file exists
if [ ! -f "$input_file" ]; then
    >&2 echo "ERROR: input file does not exist"
    exit 1;
fi

# check that the mode is valid
if [ "$mode" != "overwrite_and_backup" ]; then
    >&2 echo "ERROR: invalid mode"
    >&2 echo "As of now, only the default mode of overwrite_and_backup is supported."
    >&2 echo "So, if you specify any other mode, the script will throw an error."
    exit 1;
fi

# perform the conversion
awk \
  -v IFS="\t" -v OFS="\t" -v flag=0 \
    '{
    if ($1 ~ /^#/) {
        print $0
    } else if ($1 ~ /^>/) {
        print $0
        if ($0 ~ /rev$/) {
            flag = 5
        } else {
            flag = 0
        }
    } else {
        if (flag == 0) {
            print $1 + 4, $3, substr($4, 5, 5) "N"
        } else if (flag == 5) {
            print $1 - 4 - 5, $3, "N" substr($4, 1, 5)
        }
    }
}' "$input_file" | insert_calling_script_header "$@" > "$output_file"

# perform overwrite_and_backup
if [ "$mode" == "overwrite_and_backup" ]; then
    mv "$input_file" "$input_file".bkup
    mv "$output_file" "$input_file"
fi