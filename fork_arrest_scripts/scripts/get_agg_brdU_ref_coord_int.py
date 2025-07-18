import argparse
import sys
import pandas as pd
from io import StringIO
from DNAscentTools.modBAM_tools import get_mod_counts_per_interval

if __name__ == "__main__":

    desc = """    Per given ref-coord intervals, gets BrdU avg and total number of modified plus unmodified bases across
    all reads

    Input: Pipe in three columns w headers contig, start, end in any order.
    Data separator should be spaces or tabs. More columns can be present
    and are retained.
    
    Output is to stdout and the same columns with the same data
    are output, but we now add two new columns brdUMean and numberBrdUAndT.
    Headers are not output. Uses tab for separator.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--modBAM', type=str,
                        required=True,
                        help='modBAM file')
    parser.add_argument('--thres', type=float,
                        required=True,
                        help='(b/w 0 and 1) threshold above (below) which T is BrdU (T)')

    args = parser.parse_args()

    # read piped in data
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    inpText = sys.stdin.read()

    # load data
    dfVal = pd.read_csv(StringIO(inpText), sep='\s+',
                        comment="#", index_col=None)

    # get BrdU amounts
    modData = get_mod_counts_per_interval(args.modBAM,
                                          zip(dfVal.contig, dfVal.start, dfVal.end),
                                          'T', 'T', args.thres, args.thres)

    unmodToModRatioAndTotBases = [((k[0] + k[1]) / (k[2] + k[3] + 1e-12), k[0] + k[1] + k[2] + k[3]) for k in zip(*modData)]
    # small number is added to denominator above to rescue cases where denominator is zero.

    dfVal['brdUMean'] = [1 / (1 + k[0]) for k in unmodToModRatioAndTotBases]

    # get total number of modified and unmodified bases per interval
    dfVal['numberBrdUAndT'] = [k[1] for k in unmodToModRatioAndTotBases]

    # write to output
    print(dfVal.to_csv(index=False, header=False, float_format='%.4f',
                       sep="\t", na_rep="NA"),
          end="")
