import pysam
import itertools
from clize import run
from fork_arrest_model_tools import BrdUVsPosModelSigmoidal, BrdUVsPosModel, BrdUVsPosModelSimpleLinear


def get_model_values_at_raw_data_coords_wrap(model: str, contig: str, start: int,
                                             end: int, params: str, fasta: str) -> str:
    """ 
    Returns two tab-sep columns (without headers) coordinate, probability where probability calculated from model.

    Args:
        model: name of model
        contig: contig on reference
        start: 0-based start position on reference
        end: 0-based end position on reference
        params: model parameters delimited by '_' e.g. "0.7_0.1_500_1000" = "high_low_midpoint_width" for sigmoid
        fasta: reference genome

    Returns:
        a string containing two columns separated by tabs and row(s) separated by newlines

    Examples:
        >>> get_model_values_at_raw_data_coords_wrap("sigmoidal", "dummyI", 2, 12, "0.7_0.1_5_-2",
        ...     "sample_files/sample.fa")
        '3\\t0.261\\n4\\t0.327\\n7\\t0.539\\n8\\t0.591\\n9\\t0.628'
        >>> get_model_values_at_raw_data_coords_wrap("sigmoidal", "dummyI", 2, 12, "0.14_0.14_NA_NA",
        ...     "sample_files/sample.fa")
        '3\\t0.140\\n4\\t0.140\\n7\\t0.140\\n8\\t0.140\\n9\\t0.140'
        >>> get_model_values_at_raw_data_coords_wrap("sigmoidal", "dummyI", 2, 4, "0.6_0.0_5_2",
        ...     "sample_files/sample.fa")
        '3\\t0.439'
    """

    # obtain parameters from input and adjust offset to set x = 0 as start of fork
    try:
        p_high, p_low, offset, width = [float(k) for k in params.split("_")]
        offset -= start
    except ValueError:
        entries = params.split("_")
        try:
            # if high and low are equal and model is sigmoidal, then offset and width do not matter.
            p_high, p_low = float(entries[0]), float(entries[1])
            assert(p_high == p_low)
            assert(model == 'sigmoidal')
            offset = 1
            width = 1
        except ValueError:
            raise ValueError("Inputs are not appropriate")
        except AssertionError:
            raise ValueError("p_low and p_high must be equal")

    # perform some initial checks. Also, only one kind of model possible for now.
    if end < start or (not (0 <= p_low <= p_high <= 1)) or (not (model == 'sigmoidal')):
        raise ValueError("Inputs are not appropriate")

    # get sequence
    seq = pysam.FastaFile(fasta).fetch(contig, start, end).upper()

    # set model
    brdu_vs_pos = BrdUVsPosModelSigmoidal(1, [p_low, p_high - p_low, width])

    # return data
    return "\n".join(
        map(
            lambda x: f"{x[0] + start}\t{x[1]:.3f}",
            get_model_values_at_raw_data_coords(brdu_vs_pos, seq, offset)
        )
    )


def get_model_values_at_raw_data_coords(model: BrdUVsPosModel, seq: str, offset: float = 0) -> list[tuple[int, float]]:
    """
    Return coordinate, prob brdu pairs for the given model at As and Ts in reference sequence

    Args:
        model: name of model
        seq: reference sequence
        offset: model input = position - offset

    Returns:
        list of tuples (position, probability BrdU)

    Examples:
        >>> linear_model = BrdUVsPosModelSimpleLinear(1, [1, 2])
        >>> get_model_values_at_raw_data_coords(linear_model, 'GATACCAT', -20)
        [(1, 43), (2, 45), (3, 47), (6, 53), (7, 55)]
        >>> get_model_values_at_raw_data_coords(linear_model, 'GATACCAT')
        [(1, 3), (2, 5), (3, 7), (6, 13), (7, 15)]
    """

    positions = [k[0] for k in zip(itertools.count(), list(seq)) if k[1] == 'T' or k[1] == 'A']
    prob_brdu = model.get_model([k - offset for k in positions])

    return list(zip(positions, prob_brdu))


if __name__ == '__main__':
    # to get help, run python program_name.py --help
    run(get_model_values_at_raw_data_coords_wrap)
