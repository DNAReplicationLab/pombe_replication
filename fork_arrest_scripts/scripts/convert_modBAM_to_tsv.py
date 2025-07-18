import argparse
import sys
from DNAscentTools.modBAM_tools_additional import ModBamRecordProcessor


if __name__ == "__main__":

    desc = """    Convert modBAM file to tsv format.

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is to stdout, with four columns with headers:
    read_id, forward_read_position, ref_position, mod_qual

    Sample usage:
        samtools view -h sample.mod.bam | python <programName.py> > sample.tsv
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tag', type=str, required=False,
                        help='ChEBI code or one letter code of base modification',
                        default='T')
    parser.add_argument('--base', type=str, required=False,
                        help='The DNA base that is modified',
                        default='T')
    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # initialize modBAM record processor. NOTE: although we set a threshold, it is not used anywhere.
    modBam = ModBamRecordProcessor(0.5, args.tag, allow_non_na_mode=True, base=args.base)

    # make tsv file
    print("read_id\tforward_read_position\tref_position\tmod_qual")
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            modBam.process_modbam_line(line.strip())
            for k in modBam.mod_data_to_table():
                print(k.read_id, k.fwd_seq_pos, k.ref_pos, k.mod_qual, sep="\t")
            modBam.delete_data()
