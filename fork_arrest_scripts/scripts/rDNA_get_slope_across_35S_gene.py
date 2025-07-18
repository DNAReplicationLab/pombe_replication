import sys
import pandas as pd

if __name__ == '__main__':

    # Goal is to get the slope of analogue density vs coordinate across the 35S gene for each read.
    # Input is the file nascentAllStepListRaw.bed (or a similar name) output by our rDNA_detectSummary.py script.
    #   We use only the second, fourth and fifth columns of this file, which are the position, read_id, and step,
    #   respectively, where step refers to the difference between two windows straddling the position.
    #   We ignore any other columns and any rows that are commented out with a '#'.
    # Output is to stdout and is a tab-delimited file with the following columns and no header:
    #   position, read_id, step

    # check that inputs have been piped in
    if sys.stdin.isatty():
        print('Usage: cat <rDNA_bed_file> | python rDNA_get_slope_across_35S_gene.py > <output_file>')
        sys.exit(1)

    # Read input file from stdin
    input_data = pd.read_csv(sys.stdin, sep='\t', usecols=[1, 3, 4], header=None, comment='#')

    # Rename columns
    input_data.columns = ['position', 'read_id', 'step']

    # Create new columns
    input_data['repeat_unit_index'] = input_data['position'] // 9137
    input_data['position_per_repeat'] = input_data['position'] - (input_data['repeat_unit_index'] * 9137)

    # Filter rows
    filtered_data = input_data[(input_data['position_per_repeat'] >= 3218) & (input_data['position_per_repeat'] <= 3418)]

    # Group by 'read_id' and 'repeat_unit_index' and pick one row at random
    grouped_data = filtered_data.groupby(['read_id', 'repeat_unit_index']).sample(n=1)

    # Output the result to stdout
    grouped_data[['position', 'read_id', 'step']].to_csv(sys.stdout, sep='\t', index=False, header=False)
