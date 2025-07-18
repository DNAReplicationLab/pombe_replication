#!/bin/bash

# Written using ChatGPT

# Our goal is to remove redundant files created by our pipeline such as .detect, .forkSense, and unsorted .bam files.
# We are removing these files using the rationale that the .detect and .forkSense have been converted to mod bam
# files, and that any unsorted .bam files have been sorted and saved with a new file name and indexed with that
# filename.
# We just check if there are files with the same name but with .mod.sorted.bam or .mod.sorted.bam.bai extensions or
# .mod.left.sorted.bam and .mod.right.sorted.bam files for .forkSense files.
# If they exist, we remove the original files.
# Please read on for more information.

# Function to recursively traverse directories and identify files
traverse_directory() {
    local dir="$1"
    local output_file="$2"

    # Loop through all files and directories in the given directory
    for entry in "$dir"/*; do
        if [ -d "$entry" ]; then
            # If entry is a directory, recursively traverse it
            traverse_directory "$entry" "$output_file"
        elif [ -f "$entry" ]; then
            # If entry is a file, check if it matches the criteria
            filename=$(basename -- "$entry")
            case "$filename" in
                *.detect)
                    # If file ends with .detect, add to output file
                    # if a corresponding .mod.sorted.bam file exists
                    if [ -f "$entry".mod.sorted.bam ] && [ -f "$entry".mod.sorted.bam.bai ]; then
                        echo "rm \"$entry\"" >> "$output_file"
                    fi
                    ;;
                *.forkSense)
                  # If file ends with .forkSense, add to output file
                  # if a corresponding .mod.left.sorted.bam file and a .mod.right.sorted.bam file exists
                  if [ -f "$entry".mod.left.sorted.bam ] && [ -f "$entry".mod.left.sorted.bam.bai ] && \
                     [ -f "$entry".mod.right.sorted.bam ] && [ -f "$entry".mod.right.sorted.bam.bai ]; then
                      echo "rm \"$entry\"" >> "$output_file"
                  fi
                  ;;
                *.bam)
                    # If file ends with .bam, check if corresponding .sorted.bam exists.
                    # If it does, add to output file.
                    sorted_bam="${entry%.bam}.sorted.bam"
                    if [ -f "$sorted_bam" ]; then
                        echo "rm \"$entry\"" >> "$output_file"
                        if [ -f "$entry".bai ]; then
                            echo "rm \"${entry}.bai\"" >> "$output_file"
                        fi
                    fi
                    ;;
                *)
                    # Ignore files that don't match the criteria
                    ;;
            esac
        fi
    done
}

# Check if directory argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <directory> <output_file>"
    echo "Goal: Our pipeline creates some redundant files. This script helps identify and writes the commands to "
    echo "      remove them to the output file. The user can then examine the output file and run it if it is correct."
    echo "      We've made the deletion step a manual one to ensure no important files are deleted accidentally."
    echo "For more information, please read the comments in the script."
    echo "Example: bash $0 /path/to/directory output_file.sh"
    exit 1
fi

directory="$1"
output_file="$2"

if [ -f "$output_file" ]; then
    echo "Output file \"$output_file\" already exists. It will be overwritten."
fi

echo -e \#\!/bin/bash > "$output_file"
echo "" >> "$output_file"

# Start traversal and gather redundant files
traverse_directory "$directory" "$output_file"

# extract a list of files from the output file
file_list=$(mktemp)
grep -oP 'rm "\K[^"]+' "$output_file" > "$file_list"

# shellcheck disable=SC2046
echo "Traversal complete. File size of files to be deleted: $(du $(cat "$file_list") -cksh | tail -n 1)"
echo "Commands to remove files matching the criteria have been written to \"$output_file\"."
echo "After examining the commands, please run them using bash or sbatch."
echo "Example: e.g. sbatch -p ei-medium --time=10:00:00 \"$output_file\"."
echo "You can also use ei-short or a lower time limit; please set suitably."

# Clean up temporary files
rm "$file_list"