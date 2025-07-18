import argparse
import sys
from DNAscentTools.modBAM_tools_additional import ModBamRecordProcessor

if __name__ == "__main__":

    desc = """    Output counts of modified, unmodified, and missing bases per read id from mod bam file. 

    Input: Pipe in modBAM file contents in plain text including or excluding headers.

    Output: Four columns of space-separated data are output to stdout
        First line is 'detectIndex modified_count unmodified_count missing_count'
        Every other line is a detectIndex and three numbers whose meanings are given by the first line.
        detectIndex is a string in the format readID_contig_start_end_orn or just the readID (depending on whether
        a custom XA tag is present in the modBAM file).
        contig, start, end are positions on the reference genome and orn is fwd/rev.

    Sample usage:
        samtools view sample.mod.bam | python <programName.py> --thres 0.05
            # bases with modification probabilities above 0.05 are considered modified.
            # bases with modification probabilities below 0.05 are considered unmodified.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--thres', type=float, required=False,
                        help='(default 0.5) threshold above (below) which T is regarded as modified (unmodified)',
                        default=0.5)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # initialize variables
    modbam_record = ModBamRecordProcessor(args.thres, args.tag, True, True)

    # output data
    print("detectIndex modified_count unmodified_count missing_count")
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            modbam_record.process_modbam_line(line.rstrip())
            mod_count, unmod_count, missing_count = modbam_record.count_bases()
            print(f"{modbam_record.read_id} {mod_count} {unmod_count} {missing_count}")
