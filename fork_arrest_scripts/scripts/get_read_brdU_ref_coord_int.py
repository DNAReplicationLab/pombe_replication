import argparse
import sys
import pandas as pd
from io import StringIO
from DNAscentTools.modBAM_tools_additional import get_mod_mean_per_read_per_interval

if __name__ == "__main__":

    desc = """    Given ref-coord intervals, gets BrdU avg per read per interval

    Input: Pipe in three columns w headers contig, start, end in any order.
    Data separator should be spaces or tabs. More columns can be present
    and are retained.
    
    Output is to stdout and the same columns with the same data
    are output, but we now add three new columns read_id, mean_brdU, total num bases (T or B).
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

    # get mean BrdU amounts per read
    for row in dfVal.itertuples(index=False):
        inputData = "\t".join(str(k) for k in row)
        for read_id, brdu_mean, n_bases in get_mod_mean_per_read_per_interval(
                args.modBAM,
                row.contig, row.start, row.end,
                args.thres):
            print(f"{inputData}\t{read_id}\t{brdu_mean:.4f}\t{n_bases}")
