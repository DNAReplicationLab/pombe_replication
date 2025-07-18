import argparse
import sys
import pandas as pd
from io import StringIO

if __name__ == "__main__":

    desc = """    Reject forks whose autocorrelation length is less than or equal to a given value.

    Input: Pipe in four columns w headers detectIndex, mean_brdU,
    and start, end in any order. Data separator should be spaces or tabs.
    detectIndex is any index that uniquely marks each fork,
    mean_brdU is a numeric value, start and end refer to the boundaries
    of the window. Program assumes the data table has been
    sorted such that within each detectIndex, fork direction
    is either up or down the table.
    
    Output is to stdout and the same four columns with the same data
    are output except forks (detectIndices) whose autocorrelation length is less than
    or equal to a given value are rejected.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--autoCorrThres', type=float, required=True,
                        help='Define autocorrelation length as that window size where '
                             'autocorrelation drops below this value')
    parser.add_argument('--lenThres', type=int, required=True,
                        help='Reject forks whose autocorrelation length is less than or equal to this value')

    args = parser.parse_args()

    # read piped in data
    if not sys.stdin.isatty():
        inpText = sys.stdin.read()
    else:
        raise NotImplementedError('Please pipe in data')

    # load data
    df = pd.read_csv(StringIO(inpText), sep='\s+', comment="#", index_col=None)

    # normalize brdu and set up autocorrelation filter
    df['norm_brdu'] = df.groupby('detectIndex')['mean_brdU'].transform(
        lambda x: (x - x.mean()) / (x.std() + 1e-10))
    df_pass = df.groupby('detectIndex')['norm_brdu'].apply(
       lambda x: x.autocorr(lag=args.lenThres) > args.autoCorrThres)
    agg_df = df.merge(df_pass, left_on='detectIndex', right_index=True, suffixes=('', '_pass'))

    # drop rows where pass is False and output to stdout
    print(agg_df[agg_df['norm_brdu_pass']][['detectIndex', 'mean_brdU', 'start', 'end']].to_csv(index=False, sep='\t'),
          end="")
