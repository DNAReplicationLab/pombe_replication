import sys
from validate_bed_format import FilteredBedFileIterator


if __name__ == '__main__':

    # print usage if no input is provided
    if sys.stdin.isatty():
        print('Prints valid data lines from a BED file, raising errors if even one invalid line is found')
        print('Usage: cat <file> | python print_valid_data_bed_lines.py')
        sys.exit(1)

    try:
        for line in FilteredBedFileIterator(sys.stdin, allow_float_score=True):
            print(line.strip())
    except ValueError:
        raise ValueError('Invalid BED file format')
