import argparse
import sys
import pandas as pd
from io import StringIO
from fork_arrest_model_tools import get_pause_raw_sgm

if __name__ == "__main__":

    desc = """    Obtain candidate pause site, time from probBrdU (raw detect).

    Input: Pipe in three columns w headers detectIndex, probBrdU,
    posOnRef in any order. Data separator should be spaces or tabs.
    detectIndex is any index that uniquely marks each block of data,
    probBrdU is a numeric value, posOnRef refers to genomic coordinates
    along a reference.
    
    Output is to stdout and columns are detectIndex, pauseDuration,
    pauseSite, gof, leftCutPt with headers. pauseSite, pauseDuration
    and leftCutPt are measured in bases. 
    
    Explanation of output columns: For each detectIndex, all posOnRef
    values are considered as candidate pause sites. Program outputs
    measurements like a goodness of fit, a pause duration per pause site.
    
    Explanation of the column leftCutPt: At each candidate pause site,
    the input data is cut to produce two pieces, which each are aligned
    to the reference curve. Using the x axis of the reference curve, the limits
    of x coordinates for the left and the right pieces are  
    (leftCutPt - len_L, leftCutPt) and (leftCutPt + pauseDuration, 
    leftCutPt + pauseDuration + len_R) respectively, where len_L, len_R are
    the lengths of the two pieces. The midpoint of the reference curve
    is at x = 0.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--thres', type=float,
                        required=False, default=0.5,
                        help='(default 0.5) threshold BrdU probability at this value')

    parser.add_argument('--width', type=int,
                        required=True,
                        help='width of simple sigmoid in num of bases')

    parser.add_argument('--iter', type=int,
                        required=True,
                        help='maximum number of iterations of each cut and align step')

    parser.add_argument('--low', type=float,
                        required=False, default=0,
                        help='(default 0) low level of sigmoid')

    parser.add_argument('--tol', type=float,
                        required=True,
                        help='maximum fraction of cut and align steps that can fail per fork')

    parser.add_argument('--high', type=float,
                        required=False, default=1,
                        help='(default 1) high level of sigmoid')

    args = parser.parse_args()

    # read piped in data
    if not sys.stdin.isatty():
        inpText = sys.stdin.read()
    else:
        raise NotImplementedError("Please pipe in inputs")

    # load data
    df_val = pd.read_csv(StringIO(inpText), sep='\s+',
                         comment="#", index_col=None)

    # process user inputs
    width = args.width
    p_high = args.high
    p_low = args.low

    # loop counter
    cnt = 0

    # fit model
    for name, group in df_val.groupby('detectIndex', sort=False):

        if ('_L_' in name and width > 0) or ('_R_' in name and width < 0):
            width = -width

        x0 = group['posOnRef'].min()
        L = group['posOnRef'].max() - x0

        x_vals = [(k - x0) / L for k in group['posOnRef']]
        y_vals = [1 if k >= args.thres else 0 for k in group['probBrdU']]

        dur, gof, pos, err, left_cut_pt = get_pause_raw_sgm(x_vals, y_vals, width / L, 1,
                                                            analytical=False, cnv_thres=args.tol,
                                                            p_low=p_low, p_high=p_high,
                                                            max_iter=args.iter)

        # if pause calculation failed, read is not output
        if int(err['NoPauseFitNoConv']) > 0:
            continue

        # ensure pause sites are all close to integers
        current_pause_site_list = [k * L + x0 for k in pos]
        if not all(abs(k - round(k, 2)) < 0.01 for k in current_pause_site_list):
            raise ValueError("Non integer pause sites output!")

        # make a pandas dataframe for output
        dataForPandas = {'detectIndex': [name for k in range(len(dur))],
                         'pauseDuration': [k * L for k in dur],
                         'pauseSite': [int(k) for k in current_pause_site_list],
                         'gof': gof,
                         'leftCutPt': [k * L for k in left_cut_pt]}

        dfFit = pd.DataFrame(data=dataForPandas,
                             index=range(len(dataForPandas['detectIndex'])))

        # write data to output
        print(dfFit.to_csv(index=False, header=(cnt == 0),
                           float_format='%.6f',
                           sep=" ", na_rep="NA"),
              end="")

        # increment counter
        cnt += 1
