#!/usr/bin/env python3
# Written by ST using ChatGPT
import sys
import argparse
import csv
from DNAscentTools.modBAM_tools_additional import process_fork_index


if __name__ == "__main__":
    desc = """    Process a TSV file containing a fork index column and output only the user specified columns.

        Input: Pipe in a TSV file (tab-separated file) with no column names containing a fork index column.
        
        Options: Required are the column number to process and the format string for the output columns.
                 Format string must be comma-separated.
                 Use $1, $2, $3, etc. to refer to the columns in the input file.
                 You can specify a constant string by enclosing it in quotes.

        Output is in TSV format and is to stdout.

        Sample usage:
            cat input.tsv | python split_fork_index_from_tsv.py --col 1 '$contig, $2, $3, $read_id, "*"'
            # In the example above, the first column is the fork index column and the output columns are
            # 'contig', the second column of the input, the third column of the input, the 'read_id',
            # and a column where every entry is an asterisk. Contig and read_id are derived from the fork index column.
        """

    # ensure data is piped in
    if sys.stdin.isatty():
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # parse arguments
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--col', type=int, required=True,
                        help="The column number to process (1-based index).")
    parser.add_argument('columns', type=str, help="The format string for the output columns.")

    # process arguments
    args = parser.parse_args()

    # Prepare to read the TSV file and write the processed output
    reader = csv.reader(sys.stdin, delimiter='\t')

    for row in reader:
        # Convert the row into a dictionary with column numbers as keys
        row_dict = {f'${i + 1}': row[i] for i in range(len(row))}

        # Process the fork index column (adjusting for 0-based index in Python)
        read_id, contig, start, end, strand, fork_type, fork_start, fork_end = process_fork_index(row[args.col - 1])

        # Store the processed values in a dictionary
        processed = {
            '$read_id': read_id,
            '$contig': contig,
            '$start': start,
            '$end': end,
            '$strand': strand,
            '$fork_type': fork_type,
            '$fork_start': fork_start,
            '$fork_end': fork_end
        }

        # Prepare the final output row by evaluating the column string
        output_row = []
        for col in args.columns.split(','):
            col = col.strip()
            if col in row_dict:
                output_row.append(row_dict[col])
            elif col in processed:
                output_row.append(processed[col])
            else:
                output_row.append(col.strip('"'))

        # Write the output row
        print('\t'.join(output_row))
