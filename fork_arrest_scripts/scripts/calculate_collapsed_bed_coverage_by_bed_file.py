import argparse
import sys
import operator
from validate_bed_format import FilteredBedFileIterator
from fork_arrest_tools import intersect_two_genomic_intervals
from itertools import product

if __name__ == '__main__':

    desc = """    Given two bed files A and B, collapse A into intervals by strand and get spatial coverage by B.
    
    Usage: python calculate_collapsed_bed_coverage_by_pause_file.py <file_A.bed> <file_B.bed> 
           python calculate_collapsed_bed_coverage_by_pause_file.py --normalize-by-A-length <file_A.bed> <file_B.bed>
           
    CAUTION: Neither file must have stranded self-intersections (unstranded self-intersections are fine).
             We do not check for this; it is up to the user to enforce this. 
             If your bed file has '.' strands, they are treated as '+', so you must check for stranded
             self-intersections treating '.' as '+'.

    Input: - Specify two bed files in input.
           - Both files must be of the BED6+N format.
             Which means tab-separated columns of contig, start, end, name, score, strand, and
             any number of additional columns, in that order.
           - Comments start with # and are ignored.
           - No header with column names is allowed.
           - Can use 'stdin' or '-' in place of file_A.bed or file_B.bed (but not both) to read from standard input.
           - normalize_by_A_length: If specified, normalize each interval in A to unit length before intersecting with B
           - bins: Number of bins to use for output. Default is 100.
    
    Logic: Intersect every interval in A with every interval in B.
           If intersects are found, collapse them onto a single interval by strand.
        
    Output: - Five tab-separated columns to stdout with no header or column names:
              start_intersect, end_intersect, type_of_intersect, length_of_A_interval, length_of_B_interval
            - Start and end are relative to the "arrow tail" of intervals in A.
            - type_of_intersect can be "same_strand", "opposite_strand", or "unknown".
            - length_of_A_interval is the length of the interval in A.
            - length_of_B_interval is the length of the interval in B.
            - "Arrow tail" means the start of the interval if the strand is + or the intersection involves a '.' strand,
              and the end if the strand is - where we use "start" and "end" in the bed file coordinate sense.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file_A', help='First input bed file')
    parser.add_argument('file_B', help='Second input bed file')
    parser.add_argument('--normalize-by-A-length', action='store_true',
                        help='If specified, normalize each interval in A to unit length before intersecting with B')
    args = parser.parse_args()

    # if input parameters are not valid, print the help message
    if (not args.file_A) or (not args.file_B) or (args.file_A == args.file_B == '-') or \
            (args.file_A == args.file_B == 'stdin'):
        parser.print_help()
        exit(1)

    file_A_entries = []
    file_B_entries = []

    # open both the files
    if not (args.file_A == 'stdin' or args.file_A == '-'):
        fileA = open(args.file_A, 'r')
    else:
        fileA = sys.stdin

    if not (args.file_B == 'stdin' or args.file_B == '-'):
        fileB = open(args.file_B, 'r')
    else:
        fileB = sys.stdin

    # Read entries from both files and store them.
    # Simultaneously, store the list of contigs in each file.
    contigs = set()

    for input_file, output_list in [(fileA, file_A_entries), (fileB, file_B_entries)]:
        for line in FilteredBedFileIterator(input_file, min_six_columns=True, allow_float_score=True):
            fields = line.strip().split('\t')
            # - contig, start, end, name, score, strand are added to the output list with the strand as + or -
            # - the '.' strand is converted to '+'
            output_list.append([fields[0], int(fields[1]), int(fields[2]), fields[3], float(fields[4]),
                                fields[5] if fields[5] != '.' else '+'])
            contigs.add(fields[0])
        input_file.close()

    # Sort the bed entries
    # We are not sorting by strand here as we will shortly separate by strand
    file_A_entries.sort(key=operator.itemgetter(0, 1))
    file_B_entries.sort(key=operator.itemgetter(0, 1))

    # sort the contigs
    contig_list_sorted = sorted(contigs)

    # separate file_A_entries and file_B_entries by strand
    file_A_entries_by_strand = [[k for k in file_A_entries if k[5] == '+'], [k for k in file_A_entries if k[5] == '-']]
    file_B_entries_by_strand = [[k for k in file_B_entries if k[5] == '+'], [k for k in file_B_entries if k[5] == '-']]

    # perform intersection
    for file_A, file_B in product(file_A_entries_by_strand, file_B_entries_by_strand):

        a, b = 0, 0

        while a < len(file_A) and b < len(file_B):

            entryA = file_A[a]
            entryB = file_B[b]

            if entryA[0] != entryB[0]:
                # if contigs are not the same, increment pointer of list at the earlier contig in the sorted contig list
                if contig_list_sorted.index(entryA[0]) < contig_list_sorted.index(entryB[0]):
                    a += 1
                else:
                    b += 1
                continue
            else:
                # if the contigs are the same, proceed with the intersection
                max_start = max(entryA[1], entryB[1])
                min_end = min(entryA[2], entryB[2])

                if max_start < min_end:
                    intersect_type, st, en = intersect_two_genomic_intervals(
                        entryA[1], entryA[2], entryA[5], entryB[1], entryB[2], entryB[5],
                        norm_by_int_1_length=args.normalize_by_A_length)
                    if intersect_type == 'no_intersection':
                        continue
                    else:
                        print(f"{st}\t{en}\t{intersect_type}\t{entryA[2] - entryA[1]}\t{entryB[2] - entryB[1]}")

                # Move the pointer of the list with the smaller end position
                if entryA[2] < entryB[2]:
                    a += 1
                else:
                    b += 1
