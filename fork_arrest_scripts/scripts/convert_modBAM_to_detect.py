import argparse
import sys
from DNAscentTools.modBAM_tools_additional import modBAM_record_to_detect


if __name__ == "__main__":

    desc = """    Convert modBAM file to detect format.

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is to stdout.

    Sample usage:
        samtools view -h sample.mod.bam | python <programName.py> > sample.detect

    NOTE: * If not operated in unmapped mode, the script needs a fasta reference genome.
            So, you've to specify it through the --fasta option, or
            the initial comment block in the modbam file must
            have a line '@CO:\tcomment: #Genome <filename>' with the
            reference genome in fasta format. Both the fasta file and its
            associated index called <filename>.fai must be available.
            If both the --fasta option and the comment block are present,
            the comment block takes precedence.
          * The default modification code is T. Set a different
            modification using the --tag option.
          * If the --unmapped option is set, it takes precedence over any reference genome
            specified. In this case, modification data is output along the read coordinates
            and not the reference coordinates. This is not what the detect format is designed for,
            but it can be useful for some purposes. In this case, the contig name is set to 'unmapped'
            and the orientation is set to 'unstranded', so a detect header would look like
            '>read_id unmapped 0 read_length-1 unstranded' where read_id and read_length are the
            read name and length, respectively.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tag', type=str, required=False,
                        help='ChEBI code or one letter code of base modification',
                        default='T')
    parser.add_argument('--fasta', type=str, required=False,
                        help='Fasta reference genome file name', default='')
    parser.add_argument('--unmapped', action='store_true',
                        help='Output an unmapped detect record for each read (please see script description)')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # check that both unmapped and fasta options are not set
    if args.unmapped and args.fasta:
        parser.print_help(sys.stdout)
        raise NotImplementedError("Both --unmapped and --fasta options cannot be set together.")

    # set fasta file name
    fasta_file_name = args.fasta

    # make detect file
    for line in sys.stdin:
        if line.startswith("@CO\tcomment: #"):
            comment = line.replace("@CO\tcomment: #", "").rstrip()
            if comment.startswith("Genome "):
                fasta_file_name = comment.replace("Genome ", "")
            print("#" + comment)
        elif line.startswith("@"):
            continue
        else:
            detect_string = modBAM_record_to_detect(line.rstrip(), fasta_file_name if not args.unmapped else "",
                                                    code=args.tag)
            if detect_string:
                print(detect_string)
