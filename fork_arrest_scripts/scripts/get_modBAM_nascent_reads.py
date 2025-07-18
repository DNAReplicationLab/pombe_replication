import argparse
import sys
from DNAscentTools.modBAM_tools_additional import modBAM_record_windowed_average

if __name__ == "__main__":

    desc = """    Output read ids of possible nascent reads from modBAM file identified using windowed analogue density

    Input: Pipe in modBAM file contents in plain text including or excluding headers.

    Output is a list of read ids to stdout.
    BrdU vals can be added to output (see parameter descriptions)
    Windows are over bases where modification information is available. 
    Bases where modification probabilities are missing are ignored, not treated as NAs or as zeroes.
    Modbam must follow "?" notation which means that the modification status of a base where modification
    information is missing is 'unclear' as opposed to 'unmodified'. 

    Sample usage:
        # Couple of sample commands are shown below.
        samtools view sample.mod.bam | python <programName.py> --window 300 --nascentThres 0.05
            # measure BrdU densities in non-overlapping windows of 300 thymidines each
            # mark nascent reads using the threshold 0.05 (see parameter descriptions)
        samtools view sample.mod.bam | python <programName.py> --window 0 --nascentThres 0 --showVal
            # in the mode above, all read ids are output with mean BrdU per read id
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
    parser.add_argument('--showVal', action='store_true', required=False,
                        help='(default F) show both read id and max BrdU value in TSV format',
                        default=False)

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
            read_id, windowed_values = modBAM_record_windowed_average(line.rstrip(), args.window, args.thres,
                                                                      code=args.tag)
            read_id_rev, windowed_values_rev = modBAM_record_windowed_average(line.rstrip(), args.window, args.thres,
                                                                              code=args.tag, rev=True)
            if len(windowed_values) > 0:
                max_value = max([max(windowed_values), max(windowed_values_rev)])
                if max_value >= args.nascentThres:
                    if not args.showVal:
                        print(read_id)
                    else:
                        print(f"{read_id}\t{max_value:.6f}")
