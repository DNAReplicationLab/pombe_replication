import argparse
import sys
import os
import datetime
import pandas as pd
import numpy as np
from itertools import chain, count
from validate_bed_format import bed_6_file_from_stdin_to_sorted_pandas_dataframe
from split_bed_file_by_column_value_desc import split_bed_df_into_bins_by_column_value
from fork_arrest_tools import BedEntryWindowIterator


def split_by_value(df: pd.DataFrame, column: int = 0, num_splits: int = 0,
                   is_equal_genome_size: bool = False) -> iter:
    """ Split the bed file into the required number of pieces based on the value in the given column

    Args:
        df: pandas dataframe of the bed file
        column: column number to use for value (counting from 0 as the first column) (default: 0)
        num_splits: number of splits to make based on value (default: 0)
        is_equal_genome_size: if value split is used, do we want the split pieces to have equal genome size
                              instead of equal number of rows (default: False)

    Returns:
        iterator of tuples of pandas dataframes and prefixes
        e.g. an iterator that looks like [(df1, ["VS_1"]), (df2, ["VS_2"]), (df3, ["VS_3"]), (df4, ["VS_4"])]
        where df1, df2, df3, df4 are pandas dataframes and ["VS_1"], ["VS_2"], ["VS_3"], ["VS_4"] are prefixes
    """

    if num_splits <= 1:
        return iter([])
    else:
        if column == 4 or column > 5:
            df_list = split_bed_df_into_bins_by_column_value(df, column, num_splits,
                                                             is_equal_length_instead_of_equal_rows=is_equal_genome_size)
        else:
            raise ValueError("Invalid column number")

    return zip(df_list, (["VS_" + str(k)] for k in count(1)))


def grow_and_split_by_length_avoiding_intersects(df: pd.DataFrame, align_by: str, num_bases_split: int = 0,
                                                 num_bases_grow: int = 0, prefixes: list[str] = None,
                                                 grow_which_end: str = 'both') -> iter:
    """ Grow each bed entry by the given number of bases and split by given number of bases avoiding intersections

    Args:
        df: pandas dataframe of the bed file, must have columns contig, start, end, name, score, strand, soft_start,
            soft_end. Must have zero-length intervals i.e. start == end. The soft positions are used to specify
            limits of growing the bed file.
        align_by: Must be one of 'head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref'
        num_bases_split: number of bases to split length (default: 0)
        num_bases_grow: number of bases to grow region (default: 0)
        prefixes: list of prefixes (default: None)
        grow_which_end: Must be one of '5p', '3p', 'start', 'end' or 'both' (default: 'both').
                        Only grow the 5' or 3' or start or end of the bed file.

    Returns:
        iterator of tuples of pandas dataframes and prefixes
        e.g. an input like split_by_length(df, 100, ["A"]) will return
        an iterator that looks like [(df1, ["A","LS_1"]), (df2, ["A","LS_2"]), (df3, ["A","LS_3"])]
        where df1, df2, df3 are pandas dataframes and ["A","LS_1"], ["A","LS_2"], ["A","LS_3"] are prefixes
    """

    num_ends_to_grow = 2 if grow_which_end == 'both' else 1

    if prefixes is None:
        prefixes = []
    if num_bases_grow <= 0:
        return iter([])
    if num_bases_split <= 0:
        num_bases_split = num_ends_to_grow * num_bases_grow
        # This will ensure exactly one piece per bed line
    elif num_bases_grow % num_bases_split != 0:
        # raise an error if num_bases_grow is not divisible by num_bases_split
        raise ValueError(f"num_bases_grow ({num_bases_grow}) is not divisible by num_bases_split ({num_bases_split})")

    # grow the regions by the given number of bases
    # NOTE: these may extend beyond the reference genome; a downstream function has to trim them
    df['start'] = df[['start','strand']].apply(lambda x: x.iloc[0] - num_bases_grow
        if (grow_which_end in ['both','start'] or (grow_which_end == '5p' and x[1] == '+') or
           (grow_which_end == '3p' and x[1] == '-')) else x.iloc[0], axis=1)
    df['end'] = df[['end','strand']].apply(lambda x: x.iloc[0] + num_bases_grow
        if (grow_which_end in ['both','end'] or (grow_which_end == '5p' and x[1] == '-') or
           (grow_which_end == '3p' and x[1] == '+')) else x.iloc[0], axis=1)

    # get ready to create the split pieces by creating the start and end columns
    start = np.zeros((int(num_ends_to_grow * num_bases_grow / num_bases_split), df.shape[0]))
    end = np.zeros((int(num_ends_to_grow * num_bases_grow / num_bases_split), df.shape[0]))

    # create the start and end columns
    for k in range(df.shape[0]):
        for cnt, entry in zip(count(), BedEntryWindowIterator(df['start'].iloc[k], df['end'].iloc[k],
                                                              df['strand'].iloc[k], num_bases_split,
                                                              align_by,
                                                              (df['soft_start'].iloc[k], df['soft_end'].iloc[k]))):
            start[cnt][k] = entry[0]
            end[cnt][k] = entry[1]

    # get ready to create the split pieces
    df_list = []
    for _ in range(int(num_ends_to_grow * num_bases_grow / num_bases_split)):
        df_list.append(df.drop(["soft_start", "soft_end"], axis=1))

    for k, cnt in zip(df_list, count()):
        k['start'] = start[cnt]
        k['end'] = end[cnt]

    prefix_list = [prefixes]
    # use suffixes to differentiate between the split pieces if there are more than one
    if len(df_list) > 1:
        prefix_list = [prefixes + ["LS_" + str(k + 1)] for k in range(len(df_list))]
    return zip(df_list, prefix_list)


def add_suffix(bed_file_path: str, list_suffixes: list[str]) -> str:
    """ Add suffix to the bed file path to create a new file path

    Args:
        bed_file_path: path to the bed file
        list_suffixes: list of suffixes to add to the bed file path separated by "_"

    Returns:
        new file path
    """

    if bed_file_path.endswith(".bed"):
        if len(list_suffixes) == 0:
            return bed_file_path
        return "_".join([bed_file_path[:-4]] + list_suffixes) + ".bed"
    else:
        raise ValueError("Invalid filename")


def trim_bed_using_fasta_index(df: pd.DataFrame, fai_file: str) -> pd.DataFrame:
    """ Trim the bed file to the regions that are present in the reference genome

    Args:
        df: pandas dataframe of the bed file
        fai_file: path to the fai file of the reference genome

    Returns:
        trimmed pandas dataframe
    """
    # read the fai file
    fai_df = pd.read_csv(fai_file, sep='\t', header=None)

    fai_df = fai_df.rename(columns={0: 'seq_name', 1: 'length'})
    fai_df = fai_df[['seq_name', 'length']]
    fai_df['seq_name'] = fai_df['seq_name'].str.replace('>', '')

    # merge the bed file with the fai file
    df = df.merge(fai_df, left_on='contig', right_on='seq_name', how='inner')

    # trim the bed file using the fai file
    df['length_minus_one'] = df['length'] - 1
    df['start'] = df[['start', 'length_minus_one']].min(axis=1)
    df['end'] = df[['end', 'length']].min(axis=1)

    # ensure that start and end are not negative
    df['start'] = df['start'].clip(lower=0)
    df['end'] = df['end'].clip(lower=0)

    # check that start is less than end, or throw an error
    if (df['start'] > df['end']).any():
        raise ValueError("Trimmed start is greater than trimmed end")

    return df.drop(columns=['length', 'length_minus_one', 'seq_name'])


def main() -> None:
    parser = argparse.ArgumentParser(description="Grow and split BED file by column value and/or length\n"
                                                 "Usage: cat file.bed | python "
                                                 "grow_and_split_bed_file_by_column_value_and_length.py \n"
                                                 "Output: Makes files, and prints a json object with the filenames \n"
                                                 "Please read the header of the script for more information as the "
                                                 "options are quite complex",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--output-prefix", type=str, help="Output file prefix", required=True)
    parser.add_argument("--fai-file", type=str, help="fasta index file of the reference genome",
                        required=True)
    parser.add_argument("--align-by", default='head-to-tail',
                        choices=['head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref'],
                        help='Align by arrow head or tail and proceed from there to the other end for length splits. '
                             'Only used if length splits are requested i.e. --length-split-num-bases is set to > 0. ')
    parser.add_argument("--grow-region-num-bases", type=int,
                        help="Number of bases to grow region by.",
                        required=True)
    parser.add_argument("--value-column", type=int,
                        help="Column number to use for value splits. Only used if --num-value-splits > 0.",
                        default=0)
    parser.add_argument("--num-value-splits", type=int, help="Number of splits to make based on value",
                        default=0)
    parser.add_argument("--length-split-num-bases", type=int, help="Number of bases to split length",
                        default=0)
    parser.add_argument("--is-value-split-equal-genome-instead-of-equal-rows",
                        action='store_true', help="If value split is used, do we want the split pieces to have"
                                                  "equal genome size instead of equal number of rows. "
                                                  "NOTE: Since we use zero-size bed files, now "
                                                  "this option has been deprecated.")
    parser.add_argument("--allow-self-intersections", action='store_true',
                        help="Allow self intersections in the output bed files", default=False)
    parser.add_argument("--grow-which-end", choices = ['both','5p','3p','start','end'], required=False,
                        help="Only grow the 5' or 3' or start or end of the bed file", default='both')

    args = parser.parse_args()

    r''' In some analyses in the repo, we want to take a bed file, and want to do a combination of:
    1. Splitting the bed file by a column value (e.g. score, which could be transcription level for example)
      into a given number of pieces
    2. Growing each line in a bed file by a given number of bases, avoiding intersections
    3. Splitting each line in a bed file into many pieces such that each piece has a given number of bases
    
    Then, we can count the number of pauses and measure the null hypothesis of the number of pauses for the
    entire bed file.
    
    The counting becomes complicated if there are overlaps in the bed file, whether stranded or otherwise, as
    then we have to deal with questions of where to assign each pause.
    So, to avoid this complexity, this script only accepts bed files with zero length intervals i.e. start == end,
    and no two lines in the bed file should have the same contig and start values. 
    So, we can grow the bed file in this script ensuring that no intersections happen.
    
    Due to the reasons in the paragraph above, the only permitted operations and the only permitted orders are:
    2 alone, 1 and 2 together, 2 and 3 together, or 1 and 2 and 3 together, where the numbers refer to the list
    above. We need to do these steps together as we need to ensure that the bed file is grown in a way that
    avoids intersections with any other line of the input bed file.
    
    The user can choose which operations to perform by using the command line options.
    e.g.: If you only specify --grow-region-num-bases, then only operation 2 will be performed.
    e.g.: If you specify --num-value-splits and --grow-region-num-bases, then operations 1 and then 2 will be performed.
    
    If you do not care about self-intersections, then you can use the option --allow-self-intersections.
    '''

    # ensure that data is piped in
    if sys.stdin.isatty():
        raise NotImplementedError('Please pipe in a bed file e.g. cat file.bed | python <script_name>.py')

    # get parent directory of prefix
    parent_dir = os.path.dirname(args.output_prefix)

    # if parent directory doesn't exist, make it
    if not (parent_dir == "" or os.path.exists(parent_dir)):
        os.makedirs(parent_dir, exist_ok=True)

    # get bed file as pandas dataframe
    input_df = bed_6_file_from_stdin_to_sorted_pandas_dataframe()

    # Check that the dataframe has zero length intervals
    if not (input_df['start'] == input_df['end']).all():
        raise NotImplementedError('Please use zero length intervals in the input bed file if you request to split '
                                  'by length.')

    # ensure that the strand column is either '+' or '-'
    if not all(input_df['strand'].isin(['+', '-'])):
        raise ValueError("Strand column must be either '+' or '-'")

    # rename the "score" column to 4 because users input a numerical column here
    input_df = input_df.rename(columns={"score": 4})

    # check that grow-region-num-bases is positive
    if args.grow_region_num_bases <= 0:
        raise ValueError("grow-region-num-bases must be positive")

    # issue a warning that the is_value_split_equal_genome_instead_of_equal_rows option has been deprecated
    if args.is_value_split_equal_genome_instead_of_equal_rows:
        print("Warning: The --is-value-split-equal-genome-instead-of-equal-rows option has been deprecated. "
              "Since we use zero-size bed files, this option is no longer needed.", file=sys.stderr)

    if not args.allow_self_intersections:
        # check that every contig, start combination is unique, or throw an error
        if input_df[['contig', 'start']].duplicated().any():
            raise ValueError("Bed files must have unique contig, start values")

        # get the midpoints of each region and the previous and next regions
        input_df['soft_start'] = list(input_df.groupby(by=['contig'])['start'].rolling(2).mean())
        input_df['soft_end'] = input_df['soft_start'].shift(-1)
    else:
        input_df['soft_start'] = np.nan
        input_df['soft_end'] = np.nan

    # prepare the output objects of division list and file list
    division_list = []
    file_list = []

    # split the bed file using the value column into the required number of pieces,
    # and for each piece, save it, and split it further by the length column
    # and save the pieces
    file_name_template = args.output_prefix + ".bed"

    for df_1, division_name_1 in \
            chain.from_iterable([iter([(input_df, [])]),
                                 split_by_value(input_df, args.value_column - 1, args.num_value_splits)]):

        for df_2, division_name_2 in \
                chain.from_iterable([grow_and_split_by_length_avoiding_intersects(df_1.copy(), args.align_by,
                                                                                  0,
                                                                                  args.grow_region_num_bases,
                                                                                  division_name_1, args.grow_which_end),
                                     grow_and_split_by_length_avoiding_intersects(df_1.copy(), args.align_by,
                                                                                  args.length_split_num_bases,
                                                                                  args.grow_region_num_bases,
                                                                                  division_name_1, args.grow_which_end)
                                     if args.length_split_num_bases > 0 else iter([])
                                     ]):
            # write the given bed file,
            # and store the division and file name in the output objects
            file_name = add_suffix(file_name_template, division_name_2)
            with open(file_name, 'w') as f:
                f.write("# bed file made with grow_and_split_bed_file_by_column_value_and_length.py on "
                        f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            trim_bed_using_fasta_index(df_2, args.fai_file). \
                sort_values(by=['contig', 'start', 'end']). \
                astype({'start': 'int64', 'end': 'int64'}). \
                to_csv(file_name, sep='\t', header=False, index=False, float_format='%.3f', mode='a')
            division_list.append("_".join(division_name_2) if len(division_name_2) > 0 else "NA")
            file_list.append(file_name)

    # output division list and file list in json format
    print("[")
    for i in range(len(division_list)):
        print("  {")
        print(f'    "division": "{division_list[i]}",')
        print(f'    "bed_file": "{file_list[i]}"')
        print("  }")
        if i != len(division_list) - 1:
            print(",")
    print("]")


if __name__ == "__main__":
    main()
