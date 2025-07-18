#!/bin/bash

# goal
# ----
# Count the number of leading and lagging pauses in each pause file matching a regex pattern in a directory

# usage
# -----
# bash count_leading_lagging_pauses_from_pause_file.sh <directory_path> <regex_pattern>
# <directory_path>: path to directory containing pause files
# <regex_pattern>: regex pattern to match pause files in directory

# output
# ------
# tab-separated table to stdout with three columns (and column names): file_name, n_leading_pauses, n_lagging_pauses

# Written with ChatGPT

# load python
source load_package.sh -python

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: bash $0 <directory_path> <regex_pattern>"
    echo "Example: bash $0 /path/to/directory/ *pause*"
    echo "Calculate number of leading and lagging pauses in each file matching the regex pattern in the directory."
    echo "Output to stdout a tab-separated table with three columns (and column names): file_name, n_leading_pauses, n_lagging_pauses"
    exit 1
fi

DIR="$1"
PATTERN="$2"

# Print table header
echo -e "file_name\tn_leading_pauses\tn_lagging_pauses"

# Loop through each matching file
find "$DIR" -type f -name "$PATTERN" | while read -r file; do
    # Invoke python script to convert pause file to bed with +/- for lead/lag and process its output with awk
    result=$(< "$file" python convert_pause_file_to_bed.py --LeadLagToPlusMinus | \
        awk 'BEGIN {lead=0; lag=0}
        !/^#/ && $6 == "+" {lead++}
        !/^#/ && $6 == "-" {lag++}
        END {print lead, lag}')

    # Extract filename without the path
    filename=$(basename "$file")

    # Print the results
    echo -e "$filename\t$result"
done