import argparse
from validate_bed_format import FilteredBedFileIterator


if __name__ == '__main__':

    desc = """    Perform an exact intersect of two bed files (at least six columns each) and print the matching lines.
    
    Usage: python intersect_exact_bed_files.py <file1> <file2>

    Input: Specify two bed files in input.
           If file sizes are wildly different, the smaller file should go first.
           Both files must be of the BED6+N format with a + or - strand.
           Which means tab-separated columns of contig, start, end, name, score, strand, and
           any number of additional columns, in that order.
           Comments start with # and are ignored.
           No header with column names is allowed.
    
    Logic: Match contig, start, end, strand, and name columns between the two files.
    
    Output: Print matching lines to stdout.
            Print all columns of first file followed by the score and all columns after strand of the second file.
    
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file1', help='First input bed file')
    parser.add_argument('file2', help='Second input bed file')
    args = parser.parse_args()

    # if input parameters are not valid, print the help message
    if not args.file1 or not args.file2:
        parser.print_help()
        exit(1)

    # Open the two input files for main processing
    with open(args.file1, 'r') as file1, open(args.file2, 'r') as file2:

        # Create a dictionary to store the lines from file1 with the first four columns as the key
        file1_dict = {}
        for line in FilteredBedFileIterator(file1, min_six_columns=True, no_dot_strand=True, allow_float_score=True):
            fields = line.strip().split('\t')
            key = tuple(fields[:4] + fields[5:6])
            file1_dict[key] = tuple(fields[4:])

        # Iterate through file2 and look for matching lines in file1_dict
        for line in FilteredBedFileIterator(file2, min_six_columns=True, no_dot_strand=True, allow_float_score=True):
            fields = line.strip().split('\t')
            key = tuple(fields[:4] + fields[5:6])

            # If the key exists in file1_dict, print the matching columns
            if key in file1_dict:
                print('\t'.join(key[:-1] + file1_dict[key] + tuple(fields[4:5] + fields[6:])))
