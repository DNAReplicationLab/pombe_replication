import argparse
import sys
import pandas as pd
import numpy as np
import itertools
from io import StringIO
from fork_arrest_model_tools import process_for_model_fit, \
    BrdUVsPosModelSigmoidal

if __name__ == "__main__":

    desc = """    Fit model to brdU data avg over windows.

    Input: Pipe in four columns w headers detectIndex, mean_brdU,
    and start, end in any order. Data separator should be spaces or tabs.
    detectIndex is any index that uniquely marks each block of data,
    mean_brdU is a numeric value, start and end refer to the boundaries
    of the window. Program assumes the data table has been
    sorted such that within each detectIndex, fork direction
    is down the table.
    
    Output is to stdout and the same four columns with the same data
    are output, but we now add new rows denoted by the index 'bestFitModel',
    and a new column called 'windowIndex' w integers denoting window
    indices along a theoretical fork w index of 0 indicating start
    of a fork.

    NOTE: Although window boundaries are denoted by the words 'start' and
    'end', forks do not necessarily travel in the direction of start -> end.
    A better naming scheme would have been window_boundary_1 and
    window_boundary_2 to make the names neutral, but we have kept start
    and end for historical reasons to do with the code.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--win', type=int, required=True,
                        help='window size in bases')
    parser.add_argument('--fitFile', type=str,
                        required=True,
                        help='txt file where best-fit params etc. can be reported')
    parser.add_argument('--xMax', type=int,
                        required=True,
                        help='x coord of model runs from -xMax to +xMax')
    parser.add_argument('--winTol', type=float,
                        required=True,
                        help='window sizes not within win*(1+-winTol) are rejected')

    args = parser.parse_args()

    # read piped in data
    if not sys.stdin.isatty():
        inpText = sys.stdin.read()
    else:
        raise NotImplementedError('Please pipe in data')

    # load data
    dfVal = pd.read_csv(StringIO(inpText), sep='\s+',
                        comment="#", index_col=None)

    # prepare model
    model = BrdUVsPosModelSigmoidal(args.win, [0, 0, 0])
    candParams = list(itertools.product([2 * k / 100 for k in range(11)],
                                        [k / 100 for k in range(50, 90, 2)],
                                        [k for k in np.arange(1, 5, 0.1)]))

    # fit model
    dfFit, params = process_for_model_fit(dfVal, args.win,
                                          model, candParams, [-args.xMax, args.xMax],
                                          [args.win * (1 - args.winTol), args.win * (1 + args.winTol)])

    # write params to the file
    with open(args.fitFile, 'w') as fp:
        fp.write(str(params))
        fp.write("\n")

    # for presentation purposes, we are going to assign some labels
    # 'wrong','good_short','good_long','bad','model','aggData'
    dfFit['label'] = ''

    for name, group in dfFit.groupby('detectIndex', sort=False):

        N = len(group['mean_brdU'])
        if N < 10:
            dfFit.loc[dfFit.detectIndex == name, 'label'] = 'wrong'
        elif N < 30:
            dfFit.loc[dfFit.detectIndex == name, 'label'] = 'good_short'
        else:
            dfFit.loc[dfFit.detectIndex == name, 'label'] = 'good_long'

    dfFit.loc[dfFit.detectIndex == 'bestFitModel', 'label'] = 'model'
    dfFit.loc[dfFit.detectIndex == 'aggData', 'label'] = 'aggData'

    # write data and best fit model to output
    print(dfFit.to_csv(index=False, header=False, float_format='%.10f',
                       sep=" ", na_rep="NA"),
          end="")
