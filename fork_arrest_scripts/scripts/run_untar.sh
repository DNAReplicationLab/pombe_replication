#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J untar
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL

opDir=
fileName=

# if opDir or fileName is not set, exit
if [ -z "$opDir" ] || [ -z "$fileName" ]; then
    echo "opDir or fileName is not set"
    exit 1
fi

# make opDir if it does not exist
if ! [ -d "$opDir" ]; then
    mkdir "$opDir"
fi

# check if file with name fileName exists
if ! [ -f "$fileName" ]; then
    echo "File $fileName does not exist"
    exit 1
fi

# untar file
tar -xzf "$fileName" --directory "$opDir"