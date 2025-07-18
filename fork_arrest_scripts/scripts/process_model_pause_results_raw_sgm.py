import sys
import numpy as np
import pandas as pd
from user_input_tools import receive_pandas_and_params_from_user
from DNAscentTools.modBAM_tools_additional import process_fork_index


def calculate_x_uncertainty(x: list[float], gof: list[float], confidence_level: float,
                            gof_normalized: bool = False) -> tuple[float, float]:
    """ Calculates uncertainty in x set by the given confidence levels

    Given the goodness-of-fit landscape associated with some measured variable x, calculate the
    uncertainty in x (the best-fit value of x is not calculated here). First, we assume the
    pdf of x ~ exp(-gof). Then we calculate the cdf of x by integrating the pdf, which runs from
    0 to 1. Then, we locate those values of x where cdf(x) = (1-ci)/2 and (1+ci)/2 where ci is
    the given confidence interval level. If gof_normalized is set to True, then gof is multiplied
    by len(x) before any of the steps above are performed.

    Args:
        x: x values, assumed to be sorted in ascending order.
        gof: goodness of fit values vs x
        confidence_level: b/w 0 and 1, confidence interval levels. For example, to get 99% CI in x,
            set to 0.99.
        gof_normalized: set to true if goodness of fit normalized by length of data (see docstring as well)

    Returns:
        tuple of two floats, the lower and higher ends of the confidence interval

    Examples:
        >>> y = [1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 9/9]
        >>> [k*9 for k in calculate_x_uncertainty(y, [1, 0.99, 1.02, 1.03, 0.95, 0.98, 0.97, 1.05, 0.96], 1/3, True)]
        [4.0, 7.0]
        >>> y = np.round(np.linspace(-10, 10, 1000),2)
        >>> calculate_x_uncertainty(y, 0.5 * y**2, 0.9973)
        (-2.99, 2.99)

    """

    if gof_normalized:
        gof = [len(x) * k for k in gof]

    gof_min_indx = np.nanargmin(gof)
    gof_min = gof[gof_min_indx]
    pdf_not_normalized = np.exp([-(k - gof_min) for k in gof])
    cdf = np.nancumsum(pdf_not_normalized/np.nansum(pdf_not_normalized))
    cdf_at_gof_min = cdf[gof_min_indx]
    idx1, idx2 = np.searchsorted(cdf, [cdf_at_gof_min * (1 - confidence_level),
                                       confidence_level + cdf_at_gof_min * (1 - confidence_level)])

    if idx1 < 0 or idx2 > len(x) - 1:
        return np.nan, np.nan

    return x[idx1], x[idx2]


def process_model_pause_results_raw_sgm(df: pd.DataFrame, low: float, high: float, width: float,
                                        ci_level: float = 0.99) -> pd.DataFrame:
    """ Process results from pause fitting raw (or thresholded) probability data with sigmoids.

    Args:
        df: columns detectIndex (str), pauseDuration (float), pauseSite (int), gof (float), leftCutPt (float).
            detectIndex is any index that uniquely marks each block of data, and
                must have the format 'readID_contig_...' and contig can have underscores.
                data per detect index cannot be in non-successive rows.
            pauseSite, pauseDuration are candidate pause sites and associated duration, in bases.
                please note that every candidate pause site along the fork must be present, as we assume the gof
                associated with the last candidate pause site is the gof of the fork without any cuts,
                as one cannot cut after the last candidate site. For e.g. for a fork over the sequence
                ATTAGTC, the possible cuts are AT-TAGTC and ATT-AGTC. A cut at the last thymidine is not
                possible as there are no data points left.
            gof is goodness of fit.
            leftCutPt is where the immediate left of the candidate pause site maps to reference sigmoid's coordinates.
        low: low limit of sigmoid
        high: high limit of sigmoid
        width: width of sigmoid. Positive (negative) for right-moving (left-moving) fork. If detectIndex contains
            '_L_' or '_R_', then signs are set automatically.
        ci_level: (default 0.99) b/w 0 and 1. confidence interval for uncertainty in pause site. For example, if 95% CI
            is needed, set to 0.95.

    Returns:
        five input columns are retained, but only row with lowest gof per detect selected, and several columns added.
            New columns:
                contig, start, end, paramStrLeft, paramStrRight, leftPieceXmin, leftPieceXmax, rightPieceXmin,
                rightPieceXmax, model, low, high, width, gofNoCut, pauseSite99PctCILowerUsingExp,
                pauseSite99PctCIHigherUsingExp
            All are numeric except contig, paramStrLeft, paramStrRight and model, which are strings.
            paramStrs have the format 'high_low_offset_width' where units are bases and offset is reported
            in reference genome coordinates.

    Examples:
        >>> data = {'detectIndex': ['a_c1_100_200_fwd_L_110_120'] * 5 + ['d_c2_BB_2500_4500_rev_R_3000_3500'] * 5,
        ...         'pauseDuration': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        ...         'pauseSite': [21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
        ...         'gof': [0.02, 0.04, 0.01, 0.08, 0.1, 0.12, 0.14, 0.16, 0.012, 0.198],
        ...         'leftCutPt': [500, 401, 102, 903, 704, 205, 1006, 307, 108, 9] }
        >>> df_in = pd.DataFrame(data)
        >>> df_out = process_model_pause_results_raw_sgm(df_in, 0.2, 0.8, 10, 0.5)
        >>> ",".join(str(list(df_out[k])) for k in ['detectIndex', 'pauseDuration', 'pauseSite', 'gof', 'leftCutPt'])
        "['a_c1_100_200_fwd_L_110_120', 'd_c2_BB_2500_4500_rev_R_3000_3500'],[3, 9],[23, 29],[0.01, 0.012],[102, 108]"
        >>> ",".join(str(list(df_out[k])) for k in ['contig', 'start', 'end'])
        "['c1', 'c2_BB'],[21, 26],[25, 30]"
        >>> ",".join(str(list(df_out[k])) for k in ['paramStrLeft', 'paramStrRight'])
        "['0.8_0.2_-79_-10', '0.8_0.2_-79_10'],['0.8_0.2_-82_-10', '0.8_0.2_-88_10']"
        >>> ",".join(str(list(df_out[k])) for k in ['leftPieceXmin', 'leftPieceXmax', 'rightPieceXmin',
        ...                                         'rightPieceXmax'])
        '[100, 105],[102, 108],[105, 117],[107, 118]'
        >>> ",".join(str(list(df_out[k])) for k in ['model', 'low', 'high', 'width', 'gofNoCut'])
        "['sigmoid', 'sigmoid'],[0.2, 0.2],[0.8, 0.8],[-10, 10],[0.1, 0.198]"
        >>> ",".join(str(list(df_out[k])) for k in ['pauseSite99PctCILowerUsingExp', 'pauseSite99PctCIHigherUsingExp'])
        '[22, 28],[24, 30]'
    """

    # output variables from input columns
    name_list = []
    gof_list = []
    gof_no_cut_list = []
    dur_list = []
    pause_site_list = []
    left_cut_pt_list = []

    # processed output variables
    contig_list = []
    start_list = []
    end_list = []
    param_str_left_list = []
    param_str_right_list = []
    left_piece_x_min = []
    left_piece_x_max = []
    right_piece_x_min = []
    right_piece_x_max = []
    width_list = []
    pause_site_ci_low = []
    pause_site_ci_high = []

    # loop over detect indices
    for name, group in df.groupby('detectIndex', sort=False):

        # find best-fit row
        min_idx = group['gof'].idxmin()

        # find last row
        last_idx = group['pauseSite'].idxmax()

        # select one row per detect index from input columns
        pause_pos = group['pauseSite'][min_idx]
        left_cut_pt = group['leftCutPt'][min_idx]
        duration = group['pauseDuration'][min_idx]
        right_cut_pt = left_cut_pt + duration
        gof_min = group['gof'][min_idx]
        gof_no_cut = group['gof'][last_idx]

        # get uncertainty in pause site
        pause_site_ci_low_curr, pause_site_ci_high_curr = calculate_x_uncertainty(group['pauseSite'].tolist(),
                                                                                  group['gof'].tolist(), ci_level,
                                                                                  gof_normalized=True)

        # store results
        name_list.append(name)
        dur_list.append(duration)
        pause_site_list.append(pause_pos)
        gof_list.append(gof_min)
        left_cut_pt_list.append(left_cut_pt)
        gof_no_cut_list.append(gof_no_cut)
        pause_site_ci_low.append(pause_site_ci_low_curr)
        pause_site_ci_high.append(pause_site_ci_high_curr)

        # set contig
        if isinstance(name, str):
            _, contig, _, _, _, _, _, _ = process_fork_index(name)
            contig_list.append(contig)
        else:
            raise TypeError('detectIndex must be a string')

        # get start, end
        start_pos = group['pauseSite'].min()
        end_pos = group['pauseSite'].max()
        start_list.append(start_pos)
        end_list.append(end_pos)

        # get lengths of two pieces
        len_left = pause_pos - start_pos
        len_right = end_pos - pause_pos

        # set x coordinate limits of two pieces, where x is the x-axis of the reference sigmoid
        left_piece_x_min.append(left_cut_pt - len_left)
        left_piece_x_max.append(left_cut_pt)
        right_piece_x_min.append(right_cut_pt)
        right_piece_x_max.append(right_cut_pt + len_right)

        # set width
        if ('_L_' in name and width > 0) or ('_R_' in name and width < 0):
            width = -width
        width_list.append(width)

        # set parameter strings
        # format is 'high_low_offset_width'
        offset_left = pause_pos - left_cut_pt
        offset_right = pause_pos - right_cut_pt
        param_str_left_list.append(f"{high}_{low}_{offset_left}_{width}")
        param_str_right_list.append(f"{high}_{low}_{offset_right}_{width}")

    # make a pandas dataframe
    n_entries = len(name_list)

    op_data_frame = {
        'detectIndex': name_list,
        'pauseDuration': dur_list,
        'pauseSite': pause_site_list,
        'gof': gof_list,
        'leftCutPt': left_cut_pt_list,
        'contig': contig_list,
        'start': start_list,
        'end': end_list,
        'paramStrLeft': param_str_left_list,
        'paramStrRight': param_str_right_list,
        'leftPieceXmin': left_piece_x_min,
        'leftPieceXmax': left_piece_x_max,
        'rightPieceXmin': right_piece_x_min,
        'rightPieceXmax': right_piece_x_max,
        'model': ['sigmoid' for _ in range(n_entries)],
        'low': [low for _ in range(n_entries)],
        'high': [high for _ in range(n_entries)],
        'width': width_list,
        'gofNoCut': gof_no_cut_list,
        'pauseSite99PctCILowerUsingExp': pause_site_ci_low,
        'pauseSite99PctCIHigherUsingExp': pause_site_ci_high
    }

    return pd.DataFrame(data=op_data_frame,
                        index=range(len(op_data_frame['detectIndex'])))


if __name__ == '__main__':
    # describe the program
    desc = "Pipe in columns separated by spaces or tabs and with headers according to the description below.\n"
    desc += process_model_pause_results_raw_sgm.__doc__

    # get options and inputs from user
    df_val, args = receive_pandas_and_params_from_user([('width', int, None, 'width of simple sigmoid in num of bases'),
                                                        ('low', float, 0, '(default 0) low level of sigmoid'),
                                                        ('high', float, 1, '(default 1) high level of sigmoid')
                                                        ],
                                                       [], desc, sys.stdin,
                                                       required_columns=['detectIndex', 'pauseDuration', 'pauseSite',
                                                                         'gof', 'leftCutPt'])

    # process fit results and print best-fit measurements per index
    df_processed = process_model_pause_results_raw_sgm(df_val, args.low, args.high, args.width)

    print(df_processed.to_csv(index=False, header=True,
                              float_format='%.6f',
                              sep=" ", na_rep="NA"), end="")
