import argparse
import sys
from typing import List
import pandas as pd
from DNAscentTools.modBAM_tools_additional import ModBamRecordProcessor


def sample_read_function(probability_modbam_format: List[float]) -> float:
    """
    Sample function to be applied to each read.

    Args:
        probability_modbam_format (list): list of floats between 0 and 1.

    Returns:
        float: average of the list.

    Examples:
        >>> '%.2f' % sample_read_function([0.1, 0.2, 0.3])
        '0.20'
        >>> '%.2f' % sample_read_function([0.1, 0.2, 0.3, 0.4])
        '0.25'
        >>> sample_read_function([])
        0

    """
    if len(probability_modbam_format) == 0:
        return 0
    return sum(probability_modbam_format) / len(probability_modbam_format)


def sample_window_function(windowed_data: List[List[float]]) -> float:
    """
    Sample function to be applied to windowed data

    Args:
        windowed data: list of entries, each entry is a list of floats between 0 and 1 or is an empty list.

    Returns:
        float: sum of average of each window.

    Examples:
        >>> '%.2f' % sample_window_function([[0.1, 0.2, 0.3], [0.4, 0.5]])
        '0.65'
        >>> '%.2f' % sample_window_function([[0.1, 0.2, 0.3], [0.4, 0.5], []])
        '0.65'
        >>> '%.2f' % sample_window_function([[], [], []])
        '0.00'

    """
    return sum([sum(k) / len(k) for k in windowed_data if len(k) > 0])


def sliding_window(elements, window_size, step):
    """
    Sliding window generator.

    Args:
        elements: list of elements to be iterated over.
        window_size: size of the window.
        step: step size. If step < window_size, windows will overlap.

    Yields:
        list: list of elements in the window. Can be empty.

    Examples:
        >>> list(sliding_window([1, 2, 3, 4, 5], 3, 3))
        [[1, 2, 3]]
        >>> list(sliding_window([1, 2, 3, 4, 5, 6], 3, 3))
        [[1, 2, 3], [4, 5, 6]]
        >>> list(sliding_window([1, 2, 3, 4, 5], 3, 2))
        [[1, 2, 3], [3, 4, 5]]
        >>> list(sliding_window([1, 2, 3, 4, 5], 3, 1))
        [[1, 2, 3], [2, 3, 4], [3, 4, 5]]

    """
    if len(elements) <= window_size:
        return []
    for i in range(0, len(elements) - window_size + 1, step):
        yield elements[i:i + window_size]


if __name__ == "__main__":

    desc = """    Calculate various statistics per read in given list from mod bam file after thresholding. 

    Input: * Pipe in a mod bam file with or without headers.
           * Specify a tab-separated file with many columns and with column names in the first row.
             A column called 'read_id' must be present.
             All columns are retained in the output.
           * Other inputs are optional. See the help message.

    Output: Add a few columns to the input file and print to stdout.
            The added columns are user-defined statistics.

    Sample usage: samtools view -h reads.mod.bam |\ 
                     python calculate_modBAM_per_read_stats.py --reads-tsv reads.tsv > reads_stats.tsv
                  # NOTE: \ is used to break the command into two lines for readability.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--thres', type=float, required=False,
                        help='(default 0.5) threshold above (below) which T is regarded as modified (unmodified)',
                        default=0.5)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')
    parser.add_argument('--reads-tsv', type=str, required=True,
                        help='(required) tab-separated file with reads. see help message for details.')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # initialize variables
    modbam_record = ModBamRecordProcessor(args.thres, args.tag, False, True)

    # read in reads
    reads_df = pd.read_csv(args.reads_tsv, sep="\t")

    # check if read_id column is present
    if "read_id" not in reads_df.columns:
        raise ValueError("Column 'read_id' not found in reads file.")

    # output data
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            # process line
            modbam_record.process_modbam_line(line.rstrip())

            # check if read_id of current read is in reads_df
            if modbam_record.read_id in reads_df.read_id.values:

                # sample user defined function operating on whole read
                reads_df.loc[reads_df.read_id == modbam_record.read_id, "sample_read_stat"] = \
                    sample_read_function(modbam_record.probability_modbam_format)

                # sample user defined function operating on sliding window
                reads_df.loc[reads_df.read_id == modbam_record.read_id, "sample_windowed_stat"] = \
                    sample_window_function(list(sliding_window(modbam_record.probability_modbam_format, 300, 300)))

    # print output
    reads_df.to_csv(sys.stdout, sep="\t", index=False)
