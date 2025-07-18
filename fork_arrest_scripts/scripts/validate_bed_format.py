import sys
import re
import argparse
import uuid
import pandas as pd
from io import StringIO


class FilteredBedFileIterator:
    """    Iterate over a bed file and return only valid lines, raising an error if an invalid line is found.

    Args:
        file: iterator over lines in a bed file
        min_six_columns: bool, if True, require at least six columns (default is False)
        no_dot_strand: bool, if True, require strand to be + or - (default is False)
        allow_float_score: bool, if True, allow score to be a float (default is False)
        max_three_columns: bool, if True, allow only three columns (default is False)
        require_uuid: bool, if True, require a UUID in the name field (default is False)
        exactly_six_columns: bool, if True, require exactly six columns (default is False)

    Returns:
        line from file

    Raises:
        ValueError: if an invalid line is found
    """

    def __init__(self, file, min_six_columns=False, no_dot_strand=False, allow_float_score=False,
                 max_three_columns=False, require_uuid=False, exactly_six_columns=False):
        self.file = file
        self.bed_comment_start = ('browser', 'track', '#')
        self.min_six_columns = min_six_columns
        self.no_dot_strand = no_dot_strand
        self.allow_float_score = allow_float_score
        self.in_header = True
        self.max_three_columns = max_three_columns
        self.exactly_six_columns = exactly_six_columns
        self.require_uuid = require_uuid

    def __iter__(self):
        return self

    def __next__(self):
        while True:

            bed_line = next(self.file)
            if bed_line.startswith(self.bed_comment_start):
                if self.in_header:
                    continue
                else:
                    raise ValueError('Error: invalid bed file')

            self.in_header = False

            if validate_bed_line(bed_line, self.min_six_columns, self.no_dot_strand,
                                 self.allow_float_score, self.max_three_columns, self.require_uuid,
                                 self.exactly_six_columns):
                return bed_line
            else:
                raise ValueError('Error: invalid bed file')


def validate_bed_line(bed_line: str,
                      min_six_columns: bool = False,
                      no_dot_strand: bool = False,
                      allow_float_score: bool = False,
                      max_three_columns: bool = False,
                      require_uuid: bool = False,
                      exactly_six_columns: bool = False) -> bool:
    r""" Validate a line in a BED file.

    Args:
        bed_line: line to validate
        min_six_columns: (default False) require at least 6 columns
        no_dot_strand: (default False) disallow dot as strand
        allow_float_score: (default False) allow a floating point number as score
        max_three_columns: (default False) allow only three columns
        require_uuid: (default False) require a UUID in the name field
        exactly_six_columns: (default False) require exactly six columns

    Returns:
        True if line is valid, False otherwise

    Example:
        >>> validate_bed_line("chr1\t100\t200\tname\t0\t+\t100")
        True
        >>> validate_bed_line("chr1\t100\t200\tname with space\t0\t+\t100")
        False
        >>> validate_bed_line("chr1\t100\t200\tname\t0\t+\t100\t200\t0\t2\t100,100\t0,100")
        True
        >>> validate_bed_line("chr1\t100\t200", True)
        False
        >>> validate_bed_line("chr1\t100\t200\tname\t0\t.\t100", no_dot_strand=True)
        False
        >>> validate_bed_line("chr1\t100\t200\tblah")
        True
        >>> validate_bed_line("chr1\t100\t200\t+")
        True
        >>> validate_bed_line("chr1\t100\t200\t+\t55.66")
        False
        >>> validate_bed_line("chr1\t100\t200\t+\t55.66", allow_float_score=True)
        True
        >>> validate_bed_line("chr1\t100\t200\t+\t55.66", allow_float_score=True, max_three_columns=True)
        False
        >>> validate_bed_line("chr1\t100\t200\tblah", require_uuid=True)
        False
        >>> validate_bed_line("chr1\t100\t200", require_uuid=True)
        False
        >>> validate_bed_line("chr1\t100\t200\t" + str(uuid.uuid4()), require_uuid=True)
        True
        >>> validate_bed_line("chr1\t100\t200\tname\t0\t+\t100\t200\t0\t2\t100,100\t0,100", exactly_six_columns=True)
        False
    """

    columns = bed_line.strip().split('\t')

    if (min_six_columns or no_dot_strand) and len(columns) < 6:
        return False

    if exactly_six_columns and len(columns) != 6:
        return False

    if max_three_columns and len(columns) > 3:
        return False

    if require_uuid and len(columns) < 4:
        return False

    # if allow_float_score is True, then check if the score is a float and if so, set it to 0
    if len(columns) >= 5:
        if allow_float_score:
            try:
                float(columns[4])
            except ValueError:
                return False
            columns[4] = '0'

    if len(columns) >= 3:
        if validate_first_three_columns(columns[0], columns[1], columns[2]):
            if len(columns) >= 6:
                return validate_next_three_columns(columns[3], columns[4], columns[5], no_dot_strand, require_uuid)
            elif len(columns) == 5:
                return validate_next_three_columns(columns[3], columns[4], '.', no_dot_strand, require_uuid)
            elif len(columns) == 4:
                return validate_next_three_columns(columns[3], '0', '.', no_dot_strand, require_uuid)
            else:
                return True
        else:
            return False

    return False


def validate_first_three_columns(contig: str, start: str, end: str) -> bool:
    """ Validate the first three columns of a BED line.

    Args:
        contig: contig name
        start: start position
        end: end position

    Returns:
        True if columns are valid, False otherwise
    """
    pattern1 = r'^[\w_.]{1,255}$'
    pattern2 = r'^\d{1,20}$'

    return (re.match(pattern1, contig) is not None and
            re.match(pattern2, start) is not None and
            re.match(pattern2, end) is not None and
            0 <= int(start) <= (2 ** 64 - 1) and
            0 <= int(end) <= (2 ** 64 - 1) and
            int(start) <= int(end))


def validate_next_three_columns(name: str, score: str, strand: str,
                                no_dot_strand: bool = False, require_uuid: bool = False) -> bool:
    """ Validate the next three columns of a BED line.

    Args:
        name: name of the feature
        score: score of the feature
        strand: strand of the feature
        no_dot_strand: (default False) disallow dot as strand
        require_uuid: (default False) require a UUID in the name field

    Returns:
        True if columns are valid, False otherwise
    """

    pattern1 = r'^[\x21-\x7e]{1,255}$'
    pattern2 = r'^\d{1,4}$'
    pattern3 = r'^[-+.]$'

    if require_uuid:
        try:
            uuid.UUID(name)
        except ValueError:
            return False

    if no_dot_strand:
        pattern3 = r'^[-+]$'

    return (re.match(pattern1, name) is not None and
            re.match(pattern2, score) is not None and
            re.match(pattern3, strand) is not None and
            0 <= int(score) <= 1000)


def bed_6_file_from_stdin_to_sorted_pandas_dataframe(min_number_of_columns: int = 6,
                                                     only_allow_int_score: bool = False) -> pd.DataFrame:
    """ Read a BED file from stdin and return a sorted pandas dataframe.

    Args:
        min_number_of_columns: (default 6) minimum number of columns in the BED file
        only_allow_int_score: (default False) only allow integer scores i.e. integer entries in the fifth column

    Returns:
        A sorted pandas dataframe with the BED file contents and with the first six columns renamed to the standard

    """
    # read piped in bed file
    if not sys.stdin.isatty():

        list_of_bed_lines = []
        for bed_line in FilteredBedFileIterator(sys.stdin, min_six_columns=(min_number_of_columns == 6),
                                                allow_float_score=not only_allow_int_score):
            list_of_bed_lines.append(bed_line)

        if len(list_of_bed_lines) == 0:
            raise ValueError("No valid lines in input.")

        df = pd.read_csv(StringIO("\n".join(list_of_bed_lines)), sep='\t', header=None, on_bad_lines='error')
    else:
        raise NotImplementedError("Pipe in inputs please.")

    # rename the first six columns to the standard BED column names
    df = df.rename(columns={0: 'contig', 1: 'start', 2: 'end', 3: 'name', 4: 'score', 5: 'strand'})

    # sort
    return df.sort_values(by=['contig', 'start', 'end'])


if __name__ == '__main__':

    # print usage if no input is provided
    if sys.stdin.isatty():
        print('Usage: cat <file> | python validate_bed_format.py')
        print('Usage: cat <file> | python validate_bed_format.py --six-columns')
        print('Usage: cat <file> | python validate_bed_format.py --no-dot-strand')
        print('Usage: cat <file> | python validate_bed_format.py --allow-float-score')
        print('Can use multiple options at once if desired.')
        sys.exit(1)

    # parse command line arguments using argparse
    parser = argparse.ArgumentParser(description='Validate a BED file from stdin.'
                                                 'Usage: cat <file> | python validate_bed_format.py # with options if '
                                                 'desired')
    parser.add_argument('--six-columns', action='store_true', help='Require at least 6 columns')
    parser.add_argument('--exactly-six-columns', action='store_true', help='Require exactly 6 columns')
    parser.add_argument('--max-three-columns', action='store_true', help='Require only three columns')
    parser.add_argument('--no-dot-strand', action='store_true', help='Disallow dot \'.\' as strand')
    parser.add_argument('--allow-float-score', action='store_true', help='Allow a floating point number as score')
    parser.add_argument('--require-uuid', action='store_true', help='Require a UUID as the name of the feature')
    parser.add_argument('--single-base', action='store_true', help='Require end to be start + 1 i.e. single base')
    parser.add_argument('--zero-base', action='store_true', help='Require end to be equal to start i.e. zero base')
    parser.add_argument('--never-zero-base', action='store_true', help='Require end to be start + 1 or greater')
    parser.add_argument('--score-restrict-to-0-0pt05', action='store_true',
                        help='Require score to be between 0 and 0.05')
    args = parser.parse_args()

    # complain if incompatible options are used
    if args.single_base and args.zero_base:
        print("Error: single base and zero base options are incompatible.")
        sys.exit(1)
    if args.zero_base and args.never_zero_base:
        print("Error: zero base and never zero base options are incompatible.")
        sys.exit(1)

    if args.max_three_columns:
        if args.six_columns:
            print("Error: six columns and max three columns options are incompatible.")
            sys.exit(1)
        if args.exactly_six_columns:
            print("Error: max three columns and exactly six columns options are incompatible.")
            sys.exit(1)
        if args.require_uuid:
            print("Error: max three columns and require UUID options are incompatible.")
            sys.exit(1)

    line_count = 0

    try:
        for line in FilteredBedFileIterator(sys.stdin, args.six_columns, args.no_dot_strand,
                                            args.allow_float_score, args.max_three_columns, args.require_uuid,
                                            args.exactly_six_columns):
            split_bed_line = line.strip().split('\t')
            start = int(split_bed_line[1])
            end = int(split_bed_line[2])
            if args.single_base:
                if start + 1 != end:
                    raise ValueError
            elif args.zero_base:
                if start != end:
                    raise ValueError
            elif args.never_zero_base:
                if start == end:
                    raise ValueError
            if args.score_restrict_to_0_0pt05:
                score = float(split_bed_line[4])
                if score < 0 or score > 0.05:
                    raise ValueError
            line_count += 1
            continue
        if line_count > 0:
            print("valid")
        else:
            raise ValueError
    except ValueError:
        print("invalid")
