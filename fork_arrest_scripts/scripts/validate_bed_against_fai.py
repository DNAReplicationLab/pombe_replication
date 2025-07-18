import sys
import csv
import argparse
from fork_arrest_tools import GenomeCoverage
from validate_bed_format import FilteredBedFileIterator

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Validate a bed file against a fasta index file.')
    parser.add_argument('fai', metavar='fasta_index_file', type=str,
                        help='fasta index file location')
    parser.add_argument('--min-one-valid-line', action='store_true',
                        default=False,
                        help='return valid if at least one valid line is found')
    parser.add_argument('--check-genome-cov-to-tol', default=1_000_000_000_000, type=int,
                        help='Check that bed file covers genome within given tolerance (bp). Also see --is-stranded. '
                             'Default: 1_000_000_000_000 (1 trillion bp) i.e. so large that the check is effectively '
                             'disabled.')
    parser.add_argument('--is-stranded', action='store_true',
                        default=False,
                        help='If set, any check that involves strands will treat each strand separately. '
                             'e.g. --check-genome-cov-to-tol will check both strands separately')
    parser.add_argument('--debug', action='store_true',
                        default=False,
                        help='Print debug messages')
    parser.add_argument('--exclude-contigs', type=str,
                        default='',
                        help='(default \"\" i.e. unused) Comma-separated list of (regex-like) contigs to exclude from '
                             'check. e.g. --exclude-contigs chrM,chrX,chrY. '
                             'NOTE: If an excluded contig is found in the BED file, the check will fail. '
                             'NOTE: We say regex-like because the contig names are not treated as regexes, but '
                             '      the string is searched for the given substring i.e. chrI will match chrII as well.')
    parser.add_argument('--check-no-overlap', action='store_true',
                        default=False,
                        help='Check that no two intervals overlap (in a stranded manner if --is-stranded). ')

    # print usage if no input is provided
    if sys.stdin.isatty():
        print('Usage: cat <file> | python validate_bed_against_fai.py <fasta_index_file>')
        print('Usage: cat <file> | python validate_bed_against_fai.py --min-one-valid-line <fasta_index_file>')
        print('Replace <file> with the BED file to validate, and <fasta_index_file> with the FASTA index file.')
        print('Example: cat test.bed | python validate_bed_against_fai.py test.fai')
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # check that the genome coverage flag is not set with min-one-valid-line
    if args.check_genome_cov_to_tol < 1_000_000_000_000 and args.min_one_valid_line:
        print("error: cannot check genome coverage with min-one-valid-line")
        sys.exit(1)

    # get contigs from the fasta index file
    with open(args.fai, 'r', newline='', encoding='utf-8') as fai_file:
        contig_dict = {row[0]: int(row[1]) for row in csv.reader(fai_file, delimiter='\t')
                       if not any(k in row[0] for k in filter(None, args.exclude_contigs.split(',')))}

    # create a genome coverage object to calculate coverage
    genome_coverage = GenomeCoverage(contig_dict, args.is_stranded)

    # flag to check if the current bed line is valid
    is_current_line_valid = False

    # perform the checks
    for line in FilteredBedFileIterator(sys.stdin, allow_float_score=True, min_six_columns=args.is_stranded):
        try:
            genome_coverage.update_coverage_with_bed_line(line)
            is_current_line_valid = True
        except ValueError as e:
            is_current_line_valid = False

        if is_current_line_valid == args.min_one_valid_line:
            break

    if (is_current_line_valid and genome_coverage.is_genome_covered_within_tol(args.check_genome_cov_to_tol)
            and not (args.check_no_overlap and genome_coverage.is_overlapping_intervals())):
        print("valid")
    else:
        print("invalid")

    if args.debug:
        print(f"Covered bases: {genome_coverage.covered_bases}")
        print(f"Total bases: {genome_coverage.total_bases}")
