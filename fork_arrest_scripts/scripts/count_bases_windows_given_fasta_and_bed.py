import argparse
import sys
from validate_bed_format import FilteredBedFileIterator
from fork_arrest_tools import DNASequenceIterator


def count_bases_windows_given_fasta_and_bed(fasta_file: str, bed_file: str) -> str:
    """ Count bases in windows given a fasta file and a bed file.

    Args:
        fasta_file: path to the fasta file
        bed_file: path to the bed file

    Returns:
        string with windows in bed format and 5 additional columns of number of bases per window in order A, C, G, T, N
          where N = A + C + G + T
    """

    # set up output array
    output_bed_lines = []

    if not bed_file == 'stdin':
        fp = open(bed_file, 'r')
    else:
        fp = sys.stdin

    # iterate over the bed file
    for bed_line in FilteredBedFileIterator(fp, min_six_columns=False, allow_float_score=True):

        bed_line = bed_line.strip().split('\t')

        chrom = bed_line[0]
        start = int(bed_line[1])
        end = int(bed_line[2])

        try:
            strand = bed_line[5]
        except IndexError:
            strand = '.'

        # unknown strands are set to '+'
        strand_used = strand
        if strand == '.':
            strand_used = '+'  # default strand

        n_bases = []
        for base in ['A', 'C', 'G', 'T', 'N']:
            n_bases.append(sum(1 for _ in DNASequenceIterator(chrom, start, end, strand_used, fasta_file, base)))

        output_bed_lines.append("\t".join(bed_line + [str(k) for k in n_bases]))

    fp.close()

    # return the windows in bed format
    return "\n".join(output_bed_lines)


if __name__ == '__main__':

    # use the argparse library to receive inputs for the compulsory arguments of -i, -g
    # the -i argument is the path to the bed file
    # the -g argument is the path to the fasta file

    parser = argparse.ArgumentParser(description='Count bases (A,G,C,T,N) per window given fasta file and bed file.')
    parser.add_argument('-i', '--input-bed', help='path to the bed file. Can use stdin if piping in data',
                        required=True)
    parser.add_argument('-g', '--genome', help='path to the fasta file', required=True)
    args = parser.parse_args()

    # print the windows in bed format
    print(count_bases_windows_given_fasta_and_bed(args.genome, args.input_bed))
