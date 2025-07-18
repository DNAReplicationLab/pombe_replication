import argparse
import os
from collections.abc import Iterable
from validate_bed_format import FilteredBedFileIterator
from fork_arrest_tools import DNASequenceIterator
from itertools import cycle, count, groupby


def split_bed_entries_by_given_boundary(bed_entries):
    """ Split bed entries by a given boundary.

    Args:
        bed_entries: iterator over lines in a bed file. Must have a boundary column (column 7) that is an integer and
                     lies between the start (column 2) and end (column 3) positions.

    Returns:
        iterator over lines in a bed file where each entry has been split into two entries at each boundary position.
    """

    for bed_line in FilteredBedFileIterator(bed_entries):

        # split the bed line
        bed_fields = bed_line.strip().split('\t')

        start = int(bed_fields[1])
        end = int(bed_fields[2])

        # extract and check boundary
        try:
            boundary = int(bed_fields[6])
        except ValueError:
            raise ValueError("Invalid boundary value. Must be an integer.")

        if not (start <= boundary < end):
            raise ValueError("Boundary value must be >= start and < end.")

        # create two new bed entries and yield them
        bed_line_1 = bed_fields.copy()
        bed_line_1[2] = str(boundary)
        bed_line_2 = bed_fields.copy()
        bed_line_2[1] = str(boundary)

        yield '\t'.join(bed_line_1)
        yield '\t'.join(bed_line_2)


def make_windows_given_fasta_and_bed(fasta_file: str, bed_file_iter: Iterable[str], direction: str,
                                     num: int, base: str, allow_incomplete_windows: bool = False) -> str:
    """ Make windows given a fasta file and a bed file.

    Args:
        fasta_file: path to the fasta file
        bed_file_iter: iterator over lines in a bed file
        direction: either 'F' or 'R'
        num: number of bases to use for each window
        base: base to use for counting bases for each window, can be 'A', 'C', 'G', 'T', or 'N' (N = any base)
        allow_incomplete_windows: if True, allow windows that are shorter than num bases that arise as a result of
                                  tiling or start or end of a contig. default is False.

    Returns:
        string containing the windows in bed format
    """

    # set up output array
    output_bed_lines = []

    # check inputs
    assert (direction in ['F', 'R'])
    assert (num > 0)
    assert (base in ['A', 'C', 'G', 'T', 'N'])

    # iterate over the bed file
    for bed_line in FilteredBedFileIterator(bed_file_iter, min_six_columns=(base != 'N'), allow_float_score=True):

        bed_line = bed_line.strip().split('\t')

        chrom = bed_line[0]
        start = int(bed_line[1])
        end = int(bed_line[2])

        try:
            strand = bed_line[5]
        except IndexError:
            strand = '.'

        try:
            name = bed_line[3]
        except IndexError:
            name = ''

        try:
            score = bed_line[4]
        except IndexError:
            score = ''

        # unknown strands are set to '+' if a generic base is used,
        # otherwise, raise an error
        strand_used = strand
        if strand == '.':
            if base == 'N':
                strand_used = '+'
            else:
                raise ValueError("Have to use +/- strand in bed file for specific bases.")

        # iterate over the sequence and create windows
        # go in reverse if direction is 'R'
        dna_seq = DNASequenceIterator(chrom, start, end, strand_used, fasta_file, base)

        if direction == 'R':
            dna_seq = reversed(list(dna_seq))

        # create groupings, where each group corresponds to bases in a window
        for _, group in groupby(zip(dna_seq, count()), key=lambda x: x[1] // num):

            # get the read id, start, and end
            group_list = list(group)
            read_id = group_list[0][0][0]
            start, end = sorted([group_list[0][0][1], group_list[-1][0][1]])
            end += 1

            # check if the window is too short or too long
            if len(group_list) < num and not allow_incomplete_windows:
                continue
            elif len(group_list) > num:
                raise ValueError("Window size is larger than requested window size. Something went wrong.")

            # create the bed entry
            bed_entry = [read_id, str(start), str(end)]

            if name:
                bed_entry += [name]

            if score:
                bed_entry += [score]

            if name and score:
                bed_entry += [strand]

            output_bed_lines.append("\t".join(bed_entry))

    # return the windows in bed format
    return "\n".join(output_bed_lines)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Output windows in bed format by tiling every entry in a bed file '
                                                 'using a reference genome')
    parser.add_argument('-i', '--input-bed', help='path to the bed file', required=True)
    parser.add_argument('-g', '--genome', help='path to the fasta file', required=True)
    parser.add_argument('-n', '--num', help='number of bases to use for each window', required=True)
    parser.add_argument('-b', '--base', help='windows span given number of these, can be A, C, G, T, or N (any base)',
                        required=True)
    parser.add_argument('-d', '--direction', help='set to F, R, or FB. window boundaries are arbitrary to modulo n, '
                                                  'so they can be anchored at the start (F), end (R), or at a '
                                                  'given boundary (FB) along each bed entry. e.g.: the bed '
                                                  'entry chr_dummy 100 200 with a requested tiling of 30 bases per '
                                                  'window can be tiled as 100 130, 130 160, 160 190 (w option F) or '
                                                  '110 140, 140 170, 170 200 (w option R), or some other tiling with '
                                                  'option FB. If FB is used, the 7th column of the bed file must '
                                                  'specify a coordinate between start and end to anchor windows to.',
                        required=True)
    parser.add_argument('--allow-incomplete', help='allow incomplete windows that arise naturally due to tiling or '
                                                   'chromosome ends. '
                                                   'In the example above, a last window of 190 200 (or 100 110 if R '
                                                   'is used or some other window(s) if FB is used) would be allowed.',
                        action='store_true', default=False)
    args = parser.parse_args()

    # check if input files exist
    if not os.path.isfile(args.genome):
        raise ValueError("Genome file does not exist.")

    if not os.path.isfile(args.input_bed):
        raise ValueError("Input bed file does not exist.")

    # print the windows in bed format
    if args.direction in ['F', 'R']:
        with open(args.input_bed) as bed_file_lines:
            print(make_windows_given_fasta_and_bed(args.genome, bed_file_lines, args.direction,
                                                   int(args.num), args.base, args.allow_incomplete))
    elif args.direction == 'FB':
        with open(args.input_bed) as bed_file_lines:
            for bed_file_line, current_direction in zip(split_bed_entries_by_given_boundary(bed_file_lines),
                                                        cycle(['R', 'F'])):
                output_lines = make_windows_given_fasta_and_bed(args.genome, iter([bed_file_line]), current_direction,
                                                                int(args.num), args.base, args.allow_incomplete)
                if len(output_lines) > 0:
                    print(output_lines)
    else:
        raise ValueError("Direction must be F, R, or FB.")
