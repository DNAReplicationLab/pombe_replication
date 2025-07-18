#!/bin/bash

# goal
# -----
# execute each line in a shell script as a separate job in a slurm job array

# Usage
# -----
# sbatch --time=1:00:00 -p ei-medium --array=1-"$n" job_array_converter.sh "$job_script"
# You can change the slurm parameters, but the total number of jobs and the job script must be passed as arguments
# Each line of the job must be an independent command that can be executed in a shell script.
# See the sample input script below for an example.

# sample input script
# -------------------
# echo "Hello world"
# sleep 10
# bash limited_plot_read.sh <some arguments>
# # as you can see above, each line is a separate command and will be run as a separate job as part of a job array.
# # Only the first n lines will be executed, where n is the number input by the user.
# # Each job will be executed in the order it appears in the input script and will have the same resources as
# # every other job.

# Output
# ------
# No output is produced by this script.
# The output of each job is the output of the command that was run as part of the job.

# stop execution if any command fails
set -e

# load the job script
job_script=$1

# complain if the job script does not exist
if [ ! -f "$job_script" ]
then
  echo "Job script $job_script does not exist. Exiting." >&2
  exit 1
fi

# get the line number of the last line in the job script
n=$(wc -l < "$job_script")

# get the line number of the current job
i=${SLURM_ARRAY_TASK_ID:-1000000}

# if the current job number is the default value, exit
if [ "$i" -eq 1000000 ]
then
  echo "SLURM_ARRAY_TASK_ID is not set. Exiting." >&2
  exit 1
fi

# if the current job number is greater than the total number of jobs, exit
if [ "$i" -gt "$n" ]
then
  echo "Job number $i is greater than the total number of jobs $n. Exiting." >&2
  exit 1
fi

# get the command to run
cmd=$(sed -n "${i}p" "$job_script")

# run the command
eval "$cmd"