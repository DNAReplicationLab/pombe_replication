#!/bin/bash

# goal
# -----
# print some information about the node used for the current job for debugging purposes

# usage
#------
# bash print_node_information.sh
# Incorporate this command as a step in a slurm job script to print information about the node used for the current job.

# outputs
# -------
# A bunch of plain text information about the node used for the current job is written to standard output.

# stop execution if any command fails
set -e

echo "[LOG] Start logging slurm node information"

# print free disk space in the /tmp directory
echo "Free disk space in /tmp:"
df -hl /tmp

# print slurm parameters
echo "Slurm parameters:"
scontrol show job "$SLURM_JOB_ID"

echo "Current date and time:"
date

echo "[LOG] End logging slurm node information"