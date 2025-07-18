import random
import argparse
from fork_arrest_model_tools import map_known_pause_to_msr_pause_sgm, \
    convert_false_fork_and_fit_data_to_fastq
from process_model_pause_results_raw_sgm import calculate_x_uncertainty


def main(n, p_low=0, p_high=1, detailed_fastq=False,
         numeric=False, max_iter=None):
    """ Measure pause statistics from many forks w pauses

    Default behavior
    ================
    Print a table w headers 
    fork_name actualPauseSite actualPauseDuration
    msrPauseSite msrPauseDuration gof
    and N rows of data where N is the number of forks requested
    If detailed_fastq is set
    =======================
    print a fastq string per fork. sample format:
    >readName statistic1=val1 statistic2=val2 ...
    ATGGTAGT
    +
    !"!!~!!P
    
    Args:
        n (int): number of forks
        p_low (float): (default 0) lowest possible T->B probability
        p_high (float): (default 1) highest possible T-> probability
        detailed_fastq (bool): (default F) T/F = prnt detailed fastQ/summary
        numeric (bool): (default F) T/F = numerical/analytical method
        max_iter (int): (default None, means inf) max number of iterations for numerical method

    Returns:
        None
    """

    fork_len = int(18_000)
    width = int(3_000)

    for k in range(n):
        pause_site = random.randint(0, fork_len)
        pause_duration = random.randint(0, fork_len)

        # long durations take us beyond the limit of measurement,
        # i.e. we can detect pause sites but not measure the duration reliably.
        # So we convert them to zero i.e. non-paused forks. These can be used
        # to measure true negative and false positive rates.
        if pause_duration >= fork_len - pause_site:
            pause_duration = 0

        fork_name, pos_at_min, pause_at_min, gof_min, gof, seq, err, \
            left_cut_pt_min, pos = map_known_pause_to_msr_pause_sgm(
                fork_len, pause_site, pause_duration, width, p_low, p_high,
                numeric, max_iter)

        pause_at_min = int(round(pause_at_min))
        left_cut_pt_min = int(round(left_cut_pt_min))

        # get uncertainty in pause site
        pause_site_ci_low, pause_site_ci_high = calculate_x_uncertainty(pos, gof, 0.99, gof_normalized=True)

        fastq_str = convert_false_fork_and_fit_data_to_fastq(
            fork_name, ''.join(seq), gof,
            {
                "actualPauseSite": f"{pause_site}",
                "actualPauseDuration": f"{pause_duration}",
                "msrPauseSite": f"{pos_at_min}",
                "msrPauseDuration": f"{pause_at_min}",
                "content": "gof_data",
                "pauseSite99PctCILower": f"{pause_site_ci_low}",
                "pauseSite99PctCIHigher": f"{pause_site_ci_high}"
            } | err | {
                "leftCutPt": f"{left_cut_pt_min}",
                "gofNoCut": f"{gof[-1]:.5f}"
            })

        if not detailed_fastq:
            print(fastq_str.partition('\n')[0])
        else:
            print(fastq_str, end="")


if __name__ == '__main__':
    desc = """    Calculate f where measured pause, duration = f(actual pause, duration)
    
    Needed to get statistics associated with our pause predictions.
    Basically, given an argument N, generate N right-moving forks with
    pause sites in them and associated durations chosen randomly according to 
    internal distributions. Then, using our routines, try to locate 
    the pause and measure duration, and report it.

    Output format:
    @forkName key1=val1 key2=val2 ...
    if fastq option is not set. Otherwise a detailed fastq output
    is produced, with the sequence, and a special alphanumeric string that
    represents fitness at every site using a linear scale.

    * All units are in bases. Ensure that your inputs are in bases.
    * All forkNames are fake and are generated with random read ids.
    * Fork names follow a standard format that indicates contig, length, fork
      coordinates etc.
    """

    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-n', type=int,
                        required=True,
                        help='number of forks')

    parser.add_argument('-low', type=float,
                        required=False, default=0,
                        help='(optional) lowest possible T->B probability, default 0')

    parser.add_argument('-high', type=float,
                        required=False, default=1,
                        help='(optional) highest possible T->B probability, default 1')

    parser.add_argument('-fastq', action='store_true',
                        required=False, default=False,
                        help='(optional) output a detailed fastq format')

    parser.add_argument('-numeric', action='store_true',
                        required=False, default=False,
                        help='(optional) use numerical methods instead of analytical')

    parser.add_argument('-iter', type=int,
                        required=False, default=None,
                        help='(optional) max num iterations for numerical method')

    args = parser.parse_args()

    main(args.n, args.low, args.high, args.fastq, args.numeric,
         args.iter)
