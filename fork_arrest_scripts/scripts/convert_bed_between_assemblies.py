import os
import pysam
import sys
import regex
from validate_bed_format import FilteredBedFileIterator

if __name__ == '__main__':

    desc = """    Convert coordinates of a bed file from one assembly to another.
    WARNING: This script is for simple use cases only. It is recommended that
             you use liftOver or minimap2 or similar in general. If you use
             the --allow-three-subs option, we will allow up to 3 substitutions 
             using a package regex to do approximate matches when we cannot find
             an exact match. We are not very familiar with this package, so please
             look at the results carefully to see if they make sense.
    WARNING: We assume we will find a mapping with the same contig and on the same
             strand.

    Input: A bed file on the old assembly and two fasta files, first for the old assembly
           and second for the new assembly.

    Output: A bed file on the new assembly. Start and end coordinates are converted
            to the new assembly. If the sequence in the old bed file is not found
            in the new assembly (using a simple check),
            or is found multiple times (using a simple check),
            the start and end coordinates are set to 0.
            If no match is found with simple checks, then a more complex check
            is performed if the --allow-three-subs option is used.
            We find the best match allowing for some errors.
            So if multiple matches are found, we do not output a 0 in
            this case. If no matches are found even after this, we output
            0s in the start and the end column.
            Programs like liftOver/minimap2 are more robust because they
            use more sophisticated algorithms to find the best mapping.
    """

    # get options
    if (len(sys.argv) not in [4, 5]) or (len(sys.argv) == 5 and sys.argv[4] not in ['', '--allow-three-subs']):
        print(desc)
        print("Usage: python convert_bed_between_assemblies.py <old_bed> <old_fasta> <new_fasta> [--allow-three-subs]")
        print("NOTE: [] means optional. To specify this, remove the brackets and put the option.")
        sys.exit(1)
    else:
        old_bed = sys.argv[1]
        old_fasta = sys.argv[2]
        new_fasta = sys.argv[3]
        is_allow_three_subs = not (len(sys.argv) == 4 or (len(sys.argv) > 4 and sys.argv[4] != '--allow-three-subs'))

    # check that these files exist
    if not os.path.isfile(old_bed):
        print("Error: file %s not found" % old_bed)
        sys.exit(1)
    if not os.path.isfile(old_fasta):
        print("Error: file %s not found" % old_fasta)
        sys.exit(1)
    if not os.path.isfile(new_fasta):
        print("Error: file %s not found" % new_fasta)
        sys.exit(1)

    # open the fasta files
    old_fasta_handle = pysam.FastaFile(old_fasta)
    new_fasta_handle = pysam.FastaFile(new_fasta)

    with open(old_bed) as bed_file_lines:
        for bed_line in FilteredBedFileIterator(bed_file_lines, allow_float_score=True):

            # set flags for matches
            is_match = False
            is_multiple_match = False

            # read bed data
            bed_line = bed_line.strip().split('\t', 3)
            chrom = bed_line[0]
            start = int(bed_line[1])
            end = int(bed_line[2])

            # fetch the sequence of the old bed file and the new contig, ignoring strand
            old_bed_seq = old_fasta_handle.fetch(chrom, start, end).upper()
            new_contig_seq = new_fasta_handle.fetch(chrom).upper()

            # perform simple search for the old sequence in the new contig 
            start_new = new_contig_seq.find(old_bed_seq)
            if start_new == -1:
                start_new = 0
                end_new = 0
                is_match = False
            else:
                end_new = start_new + len(old_bed_seq)
                is_match = True

            # resume the search after the first match. If there are more matches, set the coordinates to 0
            if is_match and new_contig_seq.find(old_bed_seq, start_new + 1) != -1:
                start_new = 0
                end_new = 0
                is_multiple_match = True
            else:
                is_multiple_match = False

            # try fuzzy matches if matches are not obtained so far and if the option is set
            if is_allow_three_subs and not (is_match or is_multiple_match):
                regex_fuzzy_matches = regex.search(f'({old_bed_seq}){{s<=3}}',
                                                   new_contig_seq,regex.BESTMATCH)
                if regex_fuzzy_matches:
                    start_new, end_new = regex_fuzzy_matches.span()
                    is_match = True

            # print the new bed line
            print("%s\t%d\t%d%s" % (chrom, start_new, end_new, '' if len(bed_line) == 3 else '\t' + bed_line[3]))
