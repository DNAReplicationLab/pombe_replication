import argparse
import sys
from DNAscentTools.modBAM_tools_additional import modBAM_record_windowed_average, HeaderLines


if __name__ == "__main__":

    desc = """    Output whole entries of nascent reads from modBAM identified using windowed analogue density

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is a list of nascent modBAM entries with headers to stdout.
    Windows are over bases where modification information is available. 
    Bases where modification probabilities are missing are ignored, not treated as NAs or as zeroes.
    Modbam must follow "?" notation which means that the modification status of a base where modification
    information is missing is 'unclear' as opposed to 'unmodified'. 

    Sample usage:
        samtools view sample.mod.bam | python <programName.py> --window 300 --nascentThres 0.05
            # all modBAM entries identified as nascent reads are output
        samtools view -h sample.mod.bam | python <programName.py> --window 300 --nascentThres 0.05 |\
            samtools view -Sb -o sample2.mod.bam -
            # all modBAM entries identified as nascent reads are output and piped to samtools
            # to convert to a modBAM file named sample2.mod.bam
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--window', type=int, required=False,
                        help='(default 3000) number of thymidines per window',
                        default=3000)
    parser.add_argument('--thres', type=float, required=False,
                        help='(default 0.5) threshold above (below) which T is regarded as modified (unmodified)',
                        default=0.5)
    parser.add_argument('--nascentThres', type=float, required=False,
                        help=('(default 0.05) output read id if the mean analog density of at least one window exceeds '
                              'threshold'),
                        default=0.05)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # object to store header lines
    header = HeaderLines("get_modBAM_with_nascent_reads.py", "1.0.0")

    # output data from nascent reads
    for line in sys.stdin:
        if line.startswith("@"):
            header.add_line(line.strip())
            continue
        else:
            # output header if not already done
            if not header.has_print_action_completed():
                header.add_current_program_to_header()
                header.print_lines()

            # get windowed average in forward and reverse directions
            read_id, windowed_values = modBAM_record_windowed_average(line.rstrip(), args.window, args.thres,
                                                                      code=args.tag)
            read_id_rev, windowed_values_rev = modBAM_record_windowed_average(line.rstrip(), args.window, args.thres,
                                                                              code=args.tag, rev=True)

            # output reads if at least one window exceeds threshold
            if len(windowed_values) > 0:
                max_value = max([max(windowed_values), max(windowed_values_rev)])
                if max_value >= args.nascentThres:
                    print(line.rstrip())
