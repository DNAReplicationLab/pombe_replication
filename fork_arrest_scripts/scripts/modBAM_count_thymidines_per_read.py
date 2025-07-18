import argparse
import sys
from DNAscentTools.modBAM_tools_additional import modBAM_record_windowed_average

if __name__ == "__main__":

    desc = """    Count number of thymidines (modified + unmodified) per read id in a modBAM file.
    This may not be equal to the number of thymidines on the reference genome as some may be missing.

    Input: Pipe in modBAM file contents in plain text including or excluding headers.

    Output is a tab-separated list of read ids and count of missing thymidines to stdout without headers.

    Sample usage:
        samtools view sample.mod.bam | python <programName.py>
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # output read ids of nascent reads
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            read_id, count_T = modBAM_record_windowed_average(line.rstrip(), window_size=0, threshold=0,
                                                              code=args.tag, operation_per_win=lambda x: len(x))
            print(f"{read_id}\t{count_T[0]}")
