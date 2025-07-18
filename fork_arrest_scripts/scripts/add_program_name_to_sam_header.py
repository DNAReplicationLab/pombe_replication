import sys
import argparse
from DNAscentTools.modBAM_tools_additional import HeaderLines


if __name__ == "__main__":

    desc = """    Add a header line to sam headers indicating that the program specified on the command line has run.

    Input: Pipe in header of a sam/bam/modbam file in plain text.

    Output is in plain text and is to stdout.

    Sample usage:
        samtools view -H sample.bam | python add_program_name_to_sam_header.py program_name program_version program_call
    """

    # ensure data is piped in
    if sys.stdin.isatty():
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # parse arguments
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("program_name", help="Name of program to add to header.")
    parser.add_argument("program_version", help="Version of program to add to header.")
    parser.add_argument("program_call", help="Command line call of program to add to header.")

    # process arguments
    args = parser.parse_args()

    # object to store header lines
    header = HeaderLines(args.program_name, args.program_version, args.program_call)

    # get header lines
    for line in sys.stdin:
        if line.startswith("@"):
            header.add_line(line.strip())
            continue

    # add current program and print header lines
    header.add_current_program_to_header()
    header.print_lines()
