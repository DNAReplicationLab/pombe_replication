import argparse
import pysam
import sys

# Goal: Convert a BAM file with modified base information to a FASTA file with modified bases
#       replaced by a given letter. This is like `samtools fasta` but with modified bases.
# Usage: samtools view -h <bam_filename> | python modbam2fa.py <threshold> <base>,<mod_code>,<letter> > output.fa
#        samtools view -h <bam_filename> | python modbam2fa.py <threshold> <base>,<mod_code>,<letter> |\
#          fold -w 80 > output.fa
#        (Second usage is to wrap the output at 80 characters per line)
#        <bam_filename> is the path to the BAM file
#        <threshold> is the threshold above which bases are considered as modified, between 0.004 and 1.
#                    (0.004 corresponds to approx 1/255; we cannot deal with lower values e.g. 0).
#        <base> is the nucleotide that is being modified.
#        <mod_code> is the code used in the BAM file to represent the modified base.
#        <letter> is the letter to replace the modified base with in the FASTA file.
# NOTE: We only use the basecalled sequence from the bam file here, not the aligned sequence.
#       So, there's no reverse complementing, or processing the CIGAR string etc.
# NOTE: Bases with modification data missing are regarded as unmodified.
# NOTE: If the mod BAM line uses the '.' notation, then skipped bases are regarded as unmodified.
#       To simplify our job, we limit thresholds to above 1/255 (approx 0.004), so that these bases
#       are always marked as unmodified. If a user gives a threshold as 0, then all bases have to be
#       marked as modified, and we don't want to deal with this edge case, so we limit thresholds to
#       above 0.004.
# Output: FASTQ file with modified bases replaced by the specified character
#         Sample line (where modified T bases are replaced by 'B'):
#         >read_name
#         ACGBAACCTTGG

if __name__ == "__main__":

    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Convert BAM to FASTA with modified bases replaced by a single '
                                                 'character. We use the basecalled sequence for this, not the aligned '
                                                 'sequence. E.g.: samtools view -h  <bam_filename> | '
                                                 'python modbam2fa.py <threshold> <base>,<mod_code>,<letter> > '
                                                 'output.fa | fold -w 80 > output.fa. '
                                                 'See script for more details.')
    parser.add_argument('threshold', type=float, help='Threshold for modified bases, between (approx) 0.004'
                                                      'and 1 (0.004 corresponds to approx 1/255; we cannot deal with '
                                                      'lower values e.g. 0)')
    parser.add_argument('base_mod_code_letter', type=str,
                        help='Supply three comma-separated values: base,mod_code,letter. '
                             'base is the base that is being modified (e.g. C, T etc.). '
                             'mod_code is the code used in the BAM file to represent the modified base '
                             '(e.g. m, T etc.).'
                             'letter is the letter to replace the modified base with in the FASTA file (e.g. B).')
    # Parse the arguments
    args = parser.parse_args()
    threshold = int(args.threshold * 255)
    base, mod_code, letter = args.base_mod_code_letter.split(',')

    # check if base and letter are suitable
    if base not in ['A', 'C', 'G', 'T', 'N']:
        raise ValueError("Please provide a valid base (A, C, G, T, N)")

    if not 'A' <= letter <= 'Z':
        raise ValueError("Please provide a valid letter (A-Z)")

    # mod_code can be a single letter or a number
    if not (mod_code.isdigit() or 'A' <= mod_code <= 'Z' or 'a' <= mod_code <= 'z'):
        raise ValueError("Please provide a valid mod_code (A-Z or a-z or a number)")
    if mod_code.isdigit():
        mod_code = int(mod_code)

    # check that the threshold is between 1 and 255
    if threshold < 1 or threshold > 255:
        raise ValueError("Please provide a threshold between 0.004 and 1 (0.004 corresponds to approx 1/255; do not "
                         "use 0)")

    # check that input has been piped in
    if sys.stdin.isatty():
        raise ValueError("Please pipe in the BAM file using samtools. Remember the `-h`! "
                         "e.g.: samtools view -h <bam_filename> | python modbam2fa.py <threshold> <base>,<mod_code>,"
                         "<letter> > output.fa")

    # Open the BAM file
    sam_file = pysam.AlignmentFile(sys.stdin, "r")

    # Process each read in the BAM file
    for read in sam_file.fetch():
        seq = read.get_forward_sequence()
        mod_data = read.modified_bases_forward
        print(f">{read.query_name}")
        # Create a list from the sequence to allow mutable operations
        seq_list = list(seq)
        for k in filter(lambda x: x[1] >= threshold, mod_data[(base, 0, mod_code)]):
            seq_list[k[0]] = letter
        # Join the list back into a string
        modified_seq = ''.join(seq_list)
        print(modified_seq)

    # Close the BAM file
    sam_file.close()
