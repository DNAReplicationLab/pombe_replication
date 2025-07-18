import sys
import numpy as np
import pandas as pd
from io import StringIO
from user_input_tools import receive_pandas_and_params_from_user
from get_raw_data_from_modBAM import wrap_get_read_data_from_modBAM
from fork_arrest_model_tools import tf_fit_sgm, BernoulliNot0to1
from DNAscentTools.modBAM_tools_additional import process_fork_index


def fit_sigmoid_param_str_and_get_two_gofs(mod_bam: str, detect_index: str,
                                           width: float, low: float, high: float,
                                           param_str_left: str, param_str_right: str,
                                           pause_site: float,
                                           thres: float = 0.5,
                                           max_iter: int = 20) -> tuple[str, float, float, float, float]:
    """ Fit sigmoid to data of given index from modbam file and make parameter string,
    and calculate goodness of fit for early and late pieces.

    NOTE: This function or parts of it may have already been performed in the cut-and-align procedure.
          If so, we are doing it again as a sanity check, and because we don't want to modify that piece of code.

    Args:
        mod_bam: path to mod bam file
        detect_index: fork in format of "readID_contig_start_end_orn_dirn_start_end" where orn is fwd/rev, dirn is L/R
        width: width of sigmoid
        low: low level of sigmoid
        high: high level of sigmoid
        param_str_left: best fit parameter string for left piece of the fork (before the pause),
                        format high_low_offset_width (offset is in genomic coordinates)
        param_str_right: best fit parameter string for right piece of the fork (after the pause),
                        format high_low_offset_width (offset is in genomic coordinates)
        pause_site: best-fit pause site in genomic coordinates
        thres: threshold above which T is regarded as BrdU
        max_iter: maximum number of iterations allowed

    Returns:
        parameter string, and four goodnesses of fit (more descriptions below).
        - string format "high_low_offset_width" where offset is in genome coordinates
        - goodness of fit for one sigmoid, early piece, late piece, and sum of early and late pieces.
        - the goodness of fit for one sigmoid must match the gof value without any cuts in the output of the
          cut-and-align procedure
        - the sum of the goodness of fit for the early and late pieces must match the best-fit gof value in the output
          of the cut-and-align procedure
    """

    read_id, contig, _, _, orn, dirn, start_fork, end_fork = process_fork_index(detect_index)

    data_str = wrap_get_read_data_from_modBAM(mod_bam, read_id, contig, start_fork, end_fork, True)
    df_data = pd.read_csv(StringIO(data_str), sep="\t", names=["readId", "posOnRef", "probBrdU"])

    # check that the width parameter has the right sign
    if ('_L_' in detect_index and width > 0) or ('_R_' in detect_index and width < 0):
        raise ValueError("width parameter seems to have the wrong sign")

    x0 = df_data['posOnRef'].min()

    d_type = np.float32
    x_vals = np.array([-(k - x0) / width for k in df_data['posOnRef']], d_type)
    y_vals = np.array([1 if k >= thres else 0 for k in df_data['probBrdU']], d_type)
    n = len(y_vals)

    _, _, x0_left, _ = param_str_left.split("_")
    x_vals_left = np.array([-(k - float(x0_left)) / width for k in df_data['posOnRef'] if k < pause_site], d_type)
    y_vals_left = np.array([1 if k['probBrdU'] >= thres else 0 for k in df_data.to_dict(orient="records")
                            if k['posOnRef'] < pause_site], d_type)
    log_loss_left = -1 * BernoulliNot0to1(low, high).log_prob(y_vals_left, x_vals_left)
    gof_left = float(np.mean(log_loss_left))

    _, _, x0_right, _ = param_str_right.split("_")
    x_vals_right = np.array([-(k - float(x0_right)) / width for k in df_data['posOnRef'] if k >= pause_site], d_type)
    y_vals_right = np.array([1 if k['probBrdU'] >= thres else 0 for k in df_data.to_dict(orient="records")
                             if k['posOnRef'] >= pause_site], d_type)
    log_loss_right = -1 * BernoulliNot0to1(low, high).log_prob(y_vals_right, x_vals_right)
    gof_right = float(np.mean(log_loss_right))

    gof_verify = (gof_left * np.size(log_loss_left) + gof_right * np.size(log_loss_right)) / (np.size(log_loss_left) +
                                                                                              np.size(log_loss_right))

    [beta, _, is_converged,
     _, log_loss] = [t.numpy() for t in
                     tf_fit_sgm(
                         np.ones(n, d_type).reshape((n, 1)),
                         y_vals,
                         x_vals,
                         low,
                         high,
                         max_iter)]

    if not is_converged:
        return "", -1, gof_left, gof_right, gof_verify

    offset = width * beta[0] + x0
    return f"{high}_{low}_{offset}_{width}", float(np.mean(log_loss)), gof_left, gof_right, gof_verify


def get_max_fork_probability(mod_bam_left: str, mod_bam_right: str, detect_index: str) -> float:
    """ Get maximum fork probability from modBAM files for given detect index

    Args:
        mod_bam_left: path to mod bam file containing left fork probabilities per read per position
        mod_bam_right: path to mod bam file containing right fork probabilities per read per position
        detect_index: fork in format of "readID_contig_start_end_orn_dirn_start_end" where orn is fwd/rev, dirn is L/R

    Returns:
        maximum fork probability
    """

    read_id, contig, _, _, orn, dirn, start_fork, end_fork = process_fork_index(detect_index)

    if dirn == 'L':
        mod_bam = mod_bam_left
    elif dirn == 'R':
        mod_bam = mod_bam_right
    else:
        raise ValueError("direction seems to be incorrect")

    data_str = wrap_get_read_data_from_modBAM(mod_bam, read_id, contig, start_fork, end_fork, True)
    df_data = pd.read_csv(StringIO(data_str), sep="\t", names=["readId", "posOnRef", "probBrdU"])

    return float(df_data['probBrdU'].max())


def calculate_piece_durn_fork_start(left_piece_x_min: float, left_piece_x_max: float,
                                    right_piece_x_min: float, right_piece_x_max: float,
                                    width: float) -> tuple[float, float, float]:
    """ Calculate duration of fork start event

    Args:
        left_piece_x_min: left piece xmin
        left_piece_x_max: left piece xmax
        right_piece_x_min: right piece xmin
        right_piece_x_max: right piece xmax
        width: width of sigmoid

    Returns:
        fork start relative to sigmoid inflection point, and duration of early and late pieces 
    """

    if width > 0:
        fork_start = left_piece_x_min
        early_piece_durn = left_piece_x_max - left_piece_x_min
        late_piece_durn = right_piece_x_max - right_piece_x_min
    else:
        fork_start = -right_piece_x_max
        early_piece_durn = right_piece_x_max - right_piece_x_min
        late_piece_durn = left_piece_x_max - left_piece_x_min

    return fork_start, early_piece_durn, late_piece_durn


if __name__ == '__main__':
    desc = ("Process a table of pause results using custom functions. Pipe in many columns separated by spaces "
            "or tabs. Column headers are required. Output is to stdout and separators are tabs")

    # get options and inputs from user
    df, args = receive_pandas_and_params_from_user([("modbam", str, "", "mod bam file containing analogue "
                                                                        "probabilities per read per position"),
                                                    ("modbam_left", str, "", "mod bam file containing left fork "
                                                                             "probabilities per read per position"),
                                                    ("modbam_right", str, "", "mod bam file containing right fork "
                                                                              "probabilities per read per position"),
                                                    ("keep_cols_requested", str, "",
                                                     "columns whose calculations are requested")],
                                                   [], desc, sys.stdin,
                                                   required_columns=['pauseSite', 'detectIndex'])

    # keep track of keep columns calculated here
    keep_cols_calculated = []

    # calculations of improvements in fit
    if "gof" in df.columns and "gofNoCut" in df.columns:
        df["gofImprovement"] = df["gofNoCut"] - df["gof"]

    # calculate fork length
    if "end" in df.columns and "start" in df.columns:
        df["forkLen"] = df["end"] - df["start"]

    # mark whether to keep the pause or not
    if all([x in df.columns for x in ["pauseSite99PctCIHigherUsingExp", "pauseSite99PctCILowerUsingExp",
                                      "pauseDuration"]]) and \
            any(k in args.keep_cols_requested for k in ["keep_exp_gof", "keep_exp_gof_v2"]):
        df["keep_exp_gof_v2"] = df.apply(lambda x: 
                                         (not any(np.isnan(k) for k in (
                                            x.pauseSite99PctCIHigherUsingExp,
                                            x.pauseSite99PctCILowerUsingExp, 
                                            x.pauseDuration)))
                                         and x.pauseDuration > 0
                                         and (x.pauseSite99PctCIHigherUsingExp - x.pauseSite99PctCILowerUsingExp) / 2 <= 100, axis=1)
        keep_cols_calculated.append("keep_exp_gof_v2")

    if all([x in df.columns for x in ["leftPieceXmin", "leftPieceXmax", "rightPieceXmin", "rightPieceXmax", "width"]]):
        df[["forkStart", "earlyPieceDurn", "latePieceDurn"]] = \
            df.apply(lambda x: calculate_piece_durn_fork_start(x.leftPieceXmin, x.leftPieceXmax, x.rightPieceXmin,
                                                               x.rightPieceXmax, x.width),
                     axis=1, result_type="expand")

    # obtain no-pause measurements
    if args.modbam != "" and all([x in df.columns for x in ["detectIndex", "width", "low", "high", "paramStrLeft",
                                                            "paramStrRight", "pauseSite"]]):
        df[["paramStrNoPause", "gofNoCutVerify", "gofLeft", "gofRight", "gofVerify"]] = \
            df.apply(lambda x: fit_sigmoid_param_str_and_get_two_gofs(args.modbam, x.detectIndex, x.width,
                                                                      x.low, x.high, x.paramStrLeft, x.paramStrRight,
                                                                      x.pauseSite),
                     axis=1, result_type="expand")

    # obtain fork probability measurements
    if args.modbam_left != "" and args.modbam_right != "" and "detectIndex" in df.columns:
        df["maxForkProb"] = \
            df.apply(lambda x: get_max_fork_probability(args.modbam_left, args.modbam_right, x.detectIndex),
                     axis=1)

    # the criteria used are:
    # - forkLength threshold
    # - maximum fork probability threshold
    # - early and late piece duration thresholds
    # - mandate that the best-fit model to fork should pass through y = 0.5 i.e. BrdU density = 0.5

    if "forkLen" in df.columns and "keep_fl" in args.keep_cols_requested:
        df["keep_fl"] = df.apply(lambda x: "True" if x.forkLen >= 3000 else "False", axis=1)
        keep_cols_calculated.append("keep_fl")

    if "maxForkProb" in df.columns and "keep_fp" in args.keep_cols_requested:
        df["keep_fp"] = df.apply(lambda x: "True" if x.maxForkProb >= 0.9 else "False", axis=1)
        keep_cols_calculated.append("keep_fp")

    if "earlyPieceDurn" in df.columns and "keep_ep" in args.keep_cols_requested:
        df["keep_ep"] = df.apply(lambda x: "True" if x.earlyPieceDurn >= 200 else "False", axis=1)
        keep_cols_calculated.append("keep_ep")

    if "latePieceDurn" in df.columns and "keep_lp" in args.keep_cols_requested:
        df["keep_lp"] = df.apply(lambda x: "True" if x.latePieceDurn >= 2000 else "False", axis=1)
        keep_cols_calculated.append("keep_lp")

    if all([x in df.columns for x in ["forkStart", "high", "low", "width"]]) and "keep_fs" in args.keep_cols_requested:
        df["keep_fs"] = df.apply(
            lambda x: "True" if x.forkStart < (np.log(x.high - 0.5) - np.log(0.5 - x.low)) * np.abs(
                x.width) else "False", axis=1)
        keep_cols_calculated.append("keep_fs")

    if args.keep_cols_requested != "":
        keep_cols_requested = [x.strip() for x in args.keep_cols_requested.split(",")]
        if len([x for x in keep_cols_requested if x not in keep_cols_calculated]) > 0:
            raise ValueError("Some of the requested columns are not calculated. "
                             "calculated columns are: {} ".format(",".join(keep_cols_calculated))
                             + "requested columns are: {}".format(",".join(keep_cols_requested)))

    print(df.to_csv(index=False, header=True,
                    float_format='%.6f',
                    sep="\t", na_rep="NA"), end="")
