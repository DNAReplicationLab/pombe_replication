import argparse
import sys
import pandas as pd
import numpy as np
from io import StringIO
from fork_arrest_tools import window_data_per_index


def positive_int_or_nan(x: str):
    """Converts a string to an int if it is positive, otherwise returns np.nan"""
    try:
        n = int(x)
    except ValueError:
        return np.nan
    if n <= 0:
        return np.nan
    return n


if __name__ == "__main__":

    desc = """    Calculate mean val per window within each unique index.
    If no window size is specified, mean of all data per index is reported.

    Input: Pipe in three columns w headers detectIndex, val, posOnRef
    in any order. Data separator should be spaces or tabs.
    detectIndex is any index that uniquely marks each block of data,
    val is a numeric value, posOnRef is an integer coordinate for
    position on reference. Program assumes the data table has been
    sorted such that within each detectIndex, data appears in
    coordinates arranged in an increasing order down the table.
    
    Output is to stdout and 4 columns are present: detectIndex,
    mean_val, start, end. Data separator is spaces.
    mean_val = mean of val in a window, start and end are the
    start and end of that window in reference coordinates.
    Column headers are not output.

    NOTE: 
    * If given a threshold, val is converted to 1 or 0
    depending on whether its >= thres or not.
    * Windows are non-overlapping.
    * The last window of a read is ignored if it does not conform
    to the specified size.

    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--window', type=positive_int_or_nan, required=False,
                        help='number of entries per (non-overlapping) window. defaults to avg over entire block',
                        default=np.nan)
    parser.add_argument('--thres', type=float, required=False,
                        help=('if val >= thres, then val = 1 else val = 0.'
                              ' if threshold is not set, then val is unaltered '),
                        default=np.nan)
    parser.add_argument('--notRmDupls', required=False,
                        action='store_true',
                        help=('within an index, if multiple rows w same posOnRef are '
                              'encountered, then that index is ignored unless this '
                              'option is used. Do not use the option unless you are '
                              'an experienced user. '),
                        default=False)
    parser.add_argument('--rev', required=False,
                        action='store_true',
                        help=('windows start at last data point of each index and run '
                              'leftward. do not use with --infer'),
                        default=False)
    parser.add_argument('--infer', required=False,
                        action='store_true',
                        help=('forks marked as left- or right-moving depending on whether the detectIndex '
                              'contains "_L_" or "_R_". do not use with --rev'),
                        default=False)
    parser.add_argument('--forceWinBoundaryAtPos', type=int, required=False,
                        help='force a window boundary at the supplied genomic position. not used unless given. ',
                        default=-1)
    parser.add_argument('--halfOpen', required=False,
                        action='store_true',
                        help='output half-open intervals',
                        default=False)
    parser.add_argument('--sliding', required=False,
                        action='store_true',
                        help='use sliding windows (with no gaps) instead of non-overlapping windows',
                        default=False)

    args = parser.parse_args()

    # check arguments
    if args.infer and args.rev:
        raise ValueError("Cannot use --rev and --infer simultaneously!")

    # read piped in data
    if not sys.stdin.isatty():
        inpText = sys.stdin.read()
    else:
        raise NotImplementedError('Please pipe in inputs')

    # load data
    dfVal = pd.read_csv(StringIO(inpText), sep='\s+',
                        comment="#", index_col=None)

    # window data
    if args.rev:
        dirn = 'left'
    elif args.infer:
        dirn = 'infer'
    else:
        dirn = 'right'

    dfWin = window_data_per_index(dfVal, args.thres,
                                  args.window, args.notRmDupls, dirn, args.forceWinBoundaryAtPos,
                                  args.sliding)

    if args.halfOpen:
        dfWin['end'] = dfWin['end'] + 1

    # print results.
    print(dfWin.to_csv(index=False, header=False, float_format='%.3f',
                       columns=['detectIndex', 'mean_val', 'start', 'end'], sep=" "),
          end="")
