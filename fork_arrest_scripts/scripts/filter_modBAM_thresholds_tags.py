import sys
import argparse
from DNAscentTools.modBAM_tools_additional import ModBamFilterThreshold, HeaderLines

if __name__ == "__main__":

    desc = """    Filter modBAM files to mark bases with probabilities falling between thresholds as missing
    and simultaneously remove all modifications other than the one specified.

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is in plain text and is to stdout. Convert back to .bam if need be using samtools.

    Sample usage:
        samtools view -h sample.mod.bam | python filter_modBAM_thresholds_tags.py --thresholds 0.3,0.7 --base T --mod T 
        # in the above usage, bases with T modification probabilities between 0.3 and 0.7 will be marked as missing,
        # and any modification other than on base T with mod code T will be removed.
    """

    # ensure data is piped in
    if sys.stdin.isatty():
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # set up arguments
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--thresholds", help="Comma separated list of two thresholds for marking bases as missing. "
                        "If you set them equal, then no bases will be marked as missing.",
                        required=True)
    parser.add_argument("--base", help="Base to filter.", required=True)
    parser.add_argument("--mod", help="Modification to filter.", required=True)

    args = parser.parse_args()

    # get thresholds
    thresholds = [float(k) for k in args.thresholds.split(",")]
    low_threshold = min(thresholds)
    high_threshold = max(thresholds)

    # object to store header lines
    header = HeaderLines("filter_modBAM_thresholds_tags.py", "1.0",
                         "python filter_modBAM_thresholds_tags.py --thresholds {} --base {} --mod {}".format(
                             args.thresholds, args.base, args.mod))

    # set up filter object
    modbam_filter = ModBamFilterThreshold(low_threshold, high_threshold, args.mod, args.base)

    # process input
    for line in sys.stdin:
        if line.startswith("@"):
            header.add_line(line.strip())
            continue
        else:
            # output header if not already done
            if not header.has_print_action_completed():
                # if program has been run before, then exit
                if header.get_number_of_instances_of_current_program() > 0:
                    raise NotImplementedError("Program has been run before. Exiting.")
                header.add_current_program_to_header()
                header.print_lines()

            # process line and print
            print(modbam_filter.process_modBAM_line(line.strip()))
