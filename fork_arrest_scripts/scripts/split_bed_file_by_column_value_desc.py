import os
import argparse
import json
import pandas as pd
from validate_bed_format import bed_6_file_from_stdin_to_sorted_pandas_dataframe
from itertools import pairwise


class SplitListIndexIterator:
    """    Split a given signal into N roughly equal pieces and iterate over indices of the split points.
    NOTE: when the sum of the signal is not divisible by N, there are many ways to split the signal into N pieces.
          We choose a style where every piece adds up to approximately sum(signal) // N except for the last piece
          which has the remainder of the sum(signal) // N.
    NOTE: it is up to the user to concatenate the end of the last piece which is equal to the length of the signal.

    Args:
        signal: a list of numerical values, all must be >= 0
        num_pieces: number of desired pieces after splitting (N in the description above)

    Returns:
        An iterator over the indices of the split points.

    Examples:
        >>> list(SplitListIndexIterator([1] * 10, 3))
        [0, 3, 6]
        >>> list(SplitListIndexIterator([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 3))
        [0, 8, 11]
    """

    def __init__(self, signal: list[int], num_pieces: int):

        self.signal = signal
        self.num_pieces = num_pieces
        self.signal_sum = sum(signal)

        # check that all values in signal are >= 0
        if any([x < 0 for x in signal]):
            raise ValueError('All values in signal must be >= 0.')

        self.signal_per_file = int(self.signal_sum // num_pieces)

        self.current_index = 0
        self.num_pieces_left = num_pieces
        self.signal_running_sum = -1

    def __iter__(self):
        return self

    def __next__(self):

        while True:

            if self.current_index > len(self.signal) - 1:
                raise StopIteration

            if self.signal_running_sum >= (self.num_pieces - self.num_pieces_left + 1) * self.signal_per_file:
                # if there are more pieces left, return the current index, otherwise stop the iteration
                self.num_pieces_left -= 1
                if self.num_pieces_left == 0:
                    raise StopIteration
                return self.current_index

            elif self.signal_running_sum == -1:
                # initialize the running sum and the iterator
                self.signal_running_sum = 0
                return 0
            else:
                # increment the running sum and the iterator
                self.signal_running_sum += self.signal[self.current_index]
                self.current_index += 1


def split_bed_df_into_bins_by_column_value(df: pd.DataFrame, column_num: int,
                                           num_files: int = 0, bin_boundaries: list = None,
                                           is_ascending: bool = False,
                                           is_equal_length_instead_of_equal_rows: bool = False) -> list[pd.DataFrame]:
    """    Split a bed dataframe into bins based on a numerical column.
    Option 1: split the N rows into m files with roughly equal number of rows in each file.
    To use this, set num_files to a positive integer (m in the description above).
    If is_ascending is True, the rows will be sorted in ascending order before splitting, otherwise in descending order.
    If is_equal_length_instead_of_equal_rows is True, the files will have (approximately) equal total bed entry lengths
    instead of equal number of rows.
    Option 2: split the N rows into m files based on bin boundaries. To use this, set bin_boundaries to a list of
    boundaries in ascending order. The first file will contain data in [boundary 1, boundary 2), the second file will
    contain data in [boundary 2, boundary 3), etc. The options of is_ascending and is_equal_length_instead_of_equal_rows
    are not available with this option.
    NOTE: either option 1 or option 2 must be used, but not both.

    Args:
        df: pandas dataframe with bed entries, must have columns start, end, contig, and a column with numerical
            values to split by.
        column_num: column number to split by (0-based).
        num_files: number of files to split into.
        bin_boundaries: list of bin boundaries.
        is_ascending: if True, sort in ascending order.
        is_equal_length_instead_of_equal_rows: if True, split into files with equal total bed entry lengths.

    Returns:
        A list of pandas dataframes, each dataframe corresponding to a bin.

    """
    # prepare output list
    df_bins = []

    # check for incompatible options first
    if not bin_boundaries:
        if num_files <= 0:
            raise ValueError('Either number of files or bin boundaries must be given.')
        else:
            bin_boundaries = []
    else:
        if num_files > 0:
            raise ValueError('Both number of files and bin boundaries cannot be given.')
        if is_ascending:
            raise ValueError('Ascending option cannot be given with bin boundaries.')
        if is_equal_length_instead_of_equal_rows:
            raise ValueError('equal-length-instead-of-equal-rows option cannot be given with bin boundaries.')

    if num_files > 0:

        # sort by column after dropping rows with NaN values in the column
        df = df.dropna(subset=[column_num]).sort_values(by=column_num, ascending=is_ascending)

        if not is_equal_length_instead_of_equal_rows:

            # allocate equal number of rows to each file
            split_signal = [1] * df.shape[0]

        else:

            # allocate approximately equal length per file
            split_signal = list(df['end'] - df['start'])

        # record split indices, and add an entry equal to the end of the dataframe
        split_indices = list(SplitListIndexIterator(split_signal, num_files)) + [df.shape[0]]

        for start_end in pairwise(split_indices):
            df_bins.append(df.iloc[start_end[0]:start_end[1], :])

    elif len(bin_boundaries) > 0:

        num_bins = len(bin_boundaries) - 1

        if num_bins < 1:
            raise ValueError('Insufficient bin boundaries given.')

        for j in range(num_bins):
            df_bins.append(df[(df[column_num] >= bin_boundaries[j]) & (df[column_num] < bin_boundaries[j + 1])])

    return df_bins


if __name__ == '__main__':
    desc = """    Sort a bed file by descending order of a numerical column and split into many files and sort each file
    by contig and start position. 

    Input: Pipe in a bed file. Also needed are a column number (1-based), a prefix for the output files,
           and either a number of files to split into or a set of bin boundaries, and an option of ascending
           (default is descending). 
           - If number of files is given, the bed file will be split into that many files.
             -- By default, all files will have the same number of rows except for the last one
                which will have the remainder of the rows in case the number of rows is not divisible by the number of
                files.
             -- The user can specify equal length instead of equal number of rows by using the
                --equal-length-instead-of-equal-rows option.
                In this case, sum of bed entry lengths per file will be very close to equal (not precisely
                equal because we don't want to split any bed entry into two files and because the last file will
                have the remainder of the bed entry lengths in case the total bed entry length is not divisible by
                the number of files).
           - If bin boundaries are given, first file will contain data in [boundary 1, boundary 2),
             second file will contain data in [boundary 2, boundary 3), etc.
             * boundaries should be in ascending order.
             * [) means inclusive of first boundary and exclusive of second boundary.
             * do not use --ascending option with bin boundaries.
           - If both number of files and bin boundaries are given, you'll get an error.

    Output: - A bunch of bed files with the prefix and a number appended to the end.
            - To standard output is sent a list of dictionaries with the division number and the bed file name
              in json format.
    
    Examples:
      < input.bed python split_bed_file_by_column_value_desc.py --column 5 --prefix /output/prefix --num-files 10 
      < input.bed python split_bed_file_by_column_value_desc.py --column 5 --prefix /output/prefix \
        --bin-boundaries 10,20,30,40
      < input.bed python split_bed_file_by_column_value_desc.py --column 5 --prefix /output/prefix --num-files 10 \
        --ascending
      NOTE: \ means the command continues on the next line.
        """

    # get options
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prefix', type=str, required=True, help='Prefix for output files.')
    parser.add_argument('--column', type=int, required=True,
                        help='(1-based i.e. starts from 1) Column number to sort by.')
    parser.add_argument('--num-files', type=int, required=False, default=0,
                        help='Number of files to split into.')
    parser.add_argument('--bin-boundaries', type=str, required=False, default='',
                        help='Bin boundaries separated by commas e.g. 10,20,30,40')
    parser.add_argument('--ascending', action='store_true', required=False, default=False,
                        help='Sort ascending instead of descending.')
    parser.add_argument('--equal-length-instead-of-equal-rows', action='store_true', required=False,
                        default=False,
                        help='If given, files will have (approximately) equal total bed entry lengths '
                             'instead of equal number of rows.')

    args = parser.parse_args()

    # ensure column number is either equal to 5 or larger than 6
    if not (args.column == 5 or args.column > 6):
        raise ValueError('Column number must be either 5 or larger than 6.')

    # get bed file as pandas dataframe
    input_df = bed_6_file_from_stdin_to_sorted_pandas_dataframe()

    # rename the "score" column to 4 because users input a numerical column here
    input_df = input_df.rename(columns={"score": 4})

    # get parent directory of prefix
    parent_dir = os.path.dirname(args.prefix)

    # if parent directory doesn't exist, make it
    if not (parent_dir == "" or os.path.exists(parent_dir)):
        os.makedirs(parent_dir, exist_ok=True)

    # store a list of files, and starting and ending indices for each file
    file_list = []
    start_list = []
    end_list = []

    # prepare boundary list
    if args.bin_boundaries != '':
        input_bin_boundaries = [float(x) for x in args.bin_boundaries.split(",")]
    else:
        input_bin_boundaries = []

    # store dataframes for each bin
    output_df_bins = split_bed_df_into_bins_by_column_value(input_df, args.column - 1, args.num_files,
                                                            input_bin_boundaries, args.ascending,
                                                            args.equal_length_instead_of_equal_rows)

    # create json objects and write data to files
    for i, df_bin in zip(range(len(output_df_bins)), output_df_bins):
        file_name = f'{args.prefix}_{i + 1}.bed'
        file_list.append({"division": str(i + 1), "bed_file": file_name})
        df_bin.sort_values(by=['contig', 'start', 'end']).to_csv(file_name, sep='\t', header=None, index=False)

    # print file list to stdout after converting to json
    print(json.dumps([dict(t) for t in {tuple(d.items()) for d in file_list}], indent=4))
