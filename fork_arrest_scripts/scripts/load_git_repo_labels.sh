#!/bin/bash

# script that loads git details

# load git
source load_package.sh -git

# get commit string. if uncommitted changes,
# then append _u

COMMITSTR=$(git log --pretty=format:'%h' -n 1)
DIFFS=$(git status --porcelain)

# get the date and time
TIMENOW=$(date)
