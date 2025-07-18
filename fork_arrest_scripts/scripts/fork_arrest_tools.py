import numpy as np
import pandas as pd
import pysam
import re
import copy
from typing import Tuple


def window_data_per_index(df_val, thres=np.nan,
                          win_size=np.nan, not_rm_dupls=False, dirn='right', boundary_pos=-1,
                          sliding=False):
    """ Calculate statistics per window within each unique index.

    Args:
        df_val (pandas dataframe): must have cols val, detectIndex, posOnRef. sorted so that posOnRef
            increases down the dataframe. NOTE: we do not check if the dataframe is sorted.
        thres (float): val = 1 if exceeds threshold, zero otherwise.
            default = np.nan, which means no thresholding.
        win_size (int): size of windows. default = np.nan. if default, then
            one data-spanning window per index.
        not_rm_dupls (bool): if default of False, indices with repeated posOnRef
            are ignored. Set to true to disable.
        dirn (str): 'right' (default) or 'left' or 'infer'. Start from first (last) point
            and average in a rightward (leftward) direction. If set to 'infer', dirn
            is 'right' or 'left' for detectIndices containing '_R_' and '_L_' respectively.
        boundary_pos (int): (default -1) if given, force a window boundary at this genomic position
        sliding (bool): (default False) if True, then windows are sliding windows with no gaps.
            If False, then non-overlapping windows.

    Returns:
        df_win (pandas dataframe) with cols detectIndex, mean_val, start, end
    """

    # threshold val if the user requests it
    if not np.isnan(thres):
        def thresFn(x): return 1 if x >= thres else 0

        df_val['val'] = df_val['val'].apply(thresFn)

    return_df_list = []

    if dirn == 'right':
        index0 = 0
    elif dirn == 'left':
        index0 = -1
    elif dirn == "infer":
        index0 = None
    else:
        raise ValueError('Invalid parameters!')

    # for every detect, window the data and print measurements
    for name, group in df_val.groupby('detectIndex', sort=False):

        # if duplicate posOnRef encountered, move on to the next detect
        if not not_rm_dupls and any(group.posOnRef.duplicated()):
            continue

        if sliding:

            # window data in a sliding fashion
            df_win = pd.DataFrame({'mean_val': group['val'].rolling(win_size).mean(),
                                   'start': group['posOnRef'].rolling(win_size).min(),
                                   'end': group['posOnRef'].rolling(win_size).max()
                                   })

            # drop rows where any values are NaN
            df_win.dropna(inplace=True)

            # convert start and end to int
            df_win['start'] = df_win['start'].astype(int)
            df_win['end'] = df_win['end'].astype(int)

            # add detectIndex column
            df_win['detectIndex'] = name

            # reverse the dataframe if the user requests it
            if (dirn == 'left') or (dirn == 'infer' and '_L_' in name):
                df_win = df_win[::-1].reset_index(drop=True)
            elif (dirn == 'right') or (dirn == 'infer' and '_R_' in name):
                pass
            else:
                raise ValueError('Invalid parameters!')

            # save data
            return_df_list.append(df_win)

        else:
            # set averaging direction automatically if user requests it
            if dirn == "infer":
                if '_R_' in name:
                    index0 = 0
                elif '_L_' in name:
                    index0 = -1
                else:
                    raise ValueError('Cannot determine fork direction from index!')

            # window if the user requests it
            if not np.isnan(win_size):
                curr_win_size = win_size
                # set window boundary closest to supplied position if user requests so
                if boundary_pos >= 0:
                    group['diff'] = abs(group['posOnRef'] - boundary_pos)
                    index1 = group[['diff']].idxmin()
                    if index0 == 0:
                        avg_sign = 1
                    else:
                        avg_sign = -1
                    b = group.groupby(avg_sign * (group.index - group.index[index1])
                                      // curr_win_size)
                else:
                    b = group.groupby(np.abs(group.index - group.index[index0])
                                      // curr_win_size)
            else:
                curr_win_size = len(group.index)
                b = group.groupby([1 for _ in range(curr_win_size)])

            df_win = pd.DataFrame({'mean_val': b['val'].mean(),
                                   'start': b['posOnRef'].first(),
                                   'end': b['posOnRef'].last(),
                                   'isCountCorrect': b['posOnRef'].count() == curr_win_size})

            # reject any mis-sized windows
            if not (df_win.tail(1).isCountCorrect.any()):
                df_win.drop(df_win.tail(1).index, inplace=True)

            if not (df_win.head(1).isCountCorrect.any()):
                df_win.drop(df_win.head(1).index, inplace=True)

            # save data
            df_win.drop(['isCountCorrect'], axis=1, inplace=True)
            df_win['detectIndex'] = name
            return_df_list.append(df_win)

    return pd.concat(return_df_list)


def extract_reads_pauses(df, win, thres):
    """ (deprecated function) Extracts reads with pauses in them.

    In each detectIndex, let's say the probBrdU data is b_i, 1 <= i <= N,
    where probBrdU is the probability that the thymidine at that ref
    position is a BrdU and i is the ith row of the detectIndex in the input. 
    The detectIndex is identified to have a pause if for any j,
    mean{b}_{j to j+w-1} - mean{b}_{j+w to j+2*w-1} >= threshold, where
    w = window size and threshold are specified by user.

    Args:
        df (pandas dataframe): need cols detectIndex, probBrdU
        win (int): window size
        thres (float): thres brdu frac b/w successive winds identified as a pause

    Returns:
        list of detect indices with pauses
    """

    # subtract means of neighboring windows on a rolling basis
    df['mean_diff'] = (
        df.groupby('detectIndex', sort=False)['probBrdU'].
        transform(
            lambda s: s.rolling(2 * win).apply(
                lambda x: x[0:win].mean() - x[win:2 * win].mean() >= thres, raw=False))
    )

    # return indices w pauses
    return list(pd.unique(df[df.mean_diff == 1]['detectIndex']))


def process_param_str(x: str) -> tuple[float, float, float, float]:
    """ Process model sigmoid parameter string, of the form high_low_offset_width

    Args:
        x: parameter string

    Returns:
        high, low, offset, width or np.nan, np.nan, np.nan, np.nan if the parameter string is invalid

    Examples:
        >>> process_param_str("0.5_0.46_0.7_-0.8")
        (0.5, 0.46, 0.7, -0.8)
        >>> process_param_str("-1000_2000_-4000_-0.8")
        Traceback (most recent call last):
            ...
        AssertionError
        >>> process_param_str("mast")
        (nan, nan, nan, nan)
    """
    # check if string is of the correct format
    basic_pattern = "-?\d+(\.\d+)?"

    if not (isinstance(x, str) and re.match(rf"^{basic_pattern}_{basic_pattern}_{basic_pattern}_{basic_pattern}$", x)):
        return np.nan, np.nan, np.nan, np.nan
    else:
        high, low, offset, width = [float(y) for y in x.split("_")]
        assert 1 >= high >= low >= 0
        return high, low, offset, width


def get_contig_length(contig: str, fasta_file: str) -> int:
    """ Get contig length from a fasta file

    Args:
        contig (str): contig name
        fasta_file (str): path to fasta file

    Returns:
        int: contig length
    """
    with pysam.FastaFile(fasta_file) as fasta:
        return fasta.get_reference_length(contig)


class DNASequenceIterator:
    """Iterate over a DNA sequence

    Args:
        contig (str): contig name
        start (int): start position (0-based)
        end (int): end position (0-based, excluded)
        strand (str): strand, + or -
        fasta_file (str): path to fasta file
        base (str): base to iterate over (default 'N')

    Yields:
        tuple: (contig, position, strand)

    Examples:
        >>> for x in DNASequenceIterator('dummyI', 0, 10, '+', 'sample_files/sample.fa', 'A'):
        ...     print(x)
        ('dummyI', 0, '+')
        ('dummyI', 4, '+')
        ('dummyI', 8, '+')
        >>> for x in DNASequenceIterator('dummyI', 0, 10, '-', 'sample_files/sample.fa', 'G'):
        ...     print(x)
        ('dummyI', 2, '-')
        ('dummyI', 6, '-')
        >>> for x in DNASequenceIterator('dummyI', 20, 22, '+', 'sample_files/sample.fa', 'N'):
        ...     print(x)
        ('dummyI', 20, '+')
        ('dummyI', 21, '+')
        >>> for x in DNASequenceIterator('dummyI', 20, 20, '+', 'sample_files/sample.fa', 'N'):
        ...     print(x)
        >>> for x in DNASequenceIterator('dummyI', 20, 20, '-', 'sample_files/sample.fa', 'N'):
        ...     print(x)
    """

    def __init__(self, contig, start, end, strand, fasta_file, base='N'):
        """ Initialize """
        self.contig = contig
        self.start = max(start, 0)
        self.end = min(end, get_contig_length(contig, fasta_file))
        self.strand = strand
        self.fasta_file = fasta_file
        self.base = base.upper()

        if self.start > self.end:
            raise ValueError('start must be <= end')

        self.seq = self.load_sequence().upper()

    def __iter__(self):
        """ Return iterator """
        self.current_position = self.start
        return self

    def is_within_limits(self):
        """ Check if current position is within limits """
        return self.current_position < self.end

    def move_to_next_position(self):
        """ Move to next position """
        self.current_position += 1

    def __next__(self):
        """ Return next position """
        while self.is_within_limits():
            if self.base == 'N' or self.seq[self.current_position - self.start] == self.base:
                result = (self.contig, self.current_position, self.strand)
                self.move_to_next_position()
                return result
            self.move_to_next_position()
        raise StopIteration

    def load_sequence(self):
        """ Load sequence from fasta file """
        with pysam.FastaFile(self.fasta_file) as fasta:
            sequence = fasta.fetch(self.contig, self.start, self.end).upper()
            seg = pysam.AlignedSegment()
            seg.query_sequence = sequence
            seg.flag = 16 if self.strand == '-' else 0
        forward_seq = seg.get_forward_sequence()
        if not (forward_seq is None or forward_seq == ''):
            return seg.get_forward_sequence() if self.strand == '+' else seg.get_forward_sequence()[::-1]
        else:
            return ''


class BedEntryWindowIterator:
    """Iterate over a bed entry, stepping by given window size and along given direction

    Args:
        start (int): start position (0-based), we require start < end
        end (int): end position (0-based), we require start < end
        strand (str): strand, + or -
        window_size (int): window size
        direction (str): 'head-to-tail'/'tail-to-head'/'increasing-ref'/'decreasing-ref'. Details below:
            - proceed 'head-to-tail' or 'tail-to-head' where head and tail are end and start if strand is + and
              vice versa
            - if 'increasing-ref', then just align in the direction where coordinates increase along the reference,
              and vice versa for 'decreasing-ref'
        soft_start_stop (int, int): (default = (start, end) i.e. unused) if given, then windows are created as normal,
             but if a window lies outside these coordinates, then it is set to zero size, and if a window overlaps
             with these coordinates, then only the overlapping part is retained. The parameter is specified
             as a tuple of two integers obeying the condition start <= soft_start_stop[0] <= soft_start_stop[1] <= end.
             This is useful if we want to create regular bed windows but reject some of them due to some criterion
             like overlaps with some other region; then the criterion is expressed through soft_start_stop.

    Yields:
        tuple: (start, end)

    Examples:
        >>> for x in BedEntryWindowIterator(10, 20, '+', 4, 'tail-to-head'):
        ...     print(x)
        (10, 14)
        (14, 18)
        (18, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 10, 'tail-to-head'):
        ...     print(x)
        (10, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 10, 'tail-to-head', (15, 18)):
        ...     print(x)
        (15, 18)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 10, 'tail-to-head', (15, 28)):
        ...     print(x)
        (15, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 10, 'tail-to-head', (5, 18)):
        ...     print(x)
        (10, 18)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 4, 'tail-to-head', (10, 20)):
        ...     print(x)
        (10, 14)
        (14, 18)
        (18, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 4, 'head-to-tail'):
        ...     print(x)
        (16, 20)
        (12, 16)
        (10, 12)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'head-to-tail'):
        ...     print(x)
        (10, 14)
        (14, 18)
        (18, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'tail-to-head'):
        ...     print(x)
        (16, 20)
        (12, 16)
        (10, 12)
        >>> for x in BedEntryWindowIterator(15, 17, '-', 4, 'tail-to-head'):
        ...     print(x)
        (15, 17)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 4, 'increasing-ref'):
        ...     print(x)
        (10, 14)
        (14, 18)
        (18, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '+', 4, 'decreasing-ref'):
        ...     print(x)
        (16, 20)
        (12, 16)
        (10, 12)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'increasing-ref'):
        ...     print(x)
        (10, 14)
        (14, 18)
        (18, 20)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'increasing-ref', (17, 19)):
        ...     print(x)
        (17, 17)
        (17, 18)
        (18, 19)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'decreasing-ref'):
        ...     print(x)
        (16, 20)
        (12, 16)
        (10, 12)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'decreasing-ref', (15, 15)):
        ...     print(x)
        (15, 15)
        (15, 15)
        (15, 15)
        >>> for x in BedEntryWindowIterator(10, 20, '-', 4, 'decreasing-ref', (13, 17)):
        ...     print(x)
        (16, 17)
        (13, 16)
        (13, 13)
    """

    def __init__(self, start, end, strand, window_size, direction, soft_start_stop=(np.nan, np.nan)):
        """ Initialize """

        self.start = start
        self.end = end
        self.strand = strand
        self.window_size = window_size
        self.direction = direction

        # set the soft start, soft stop coordinates
        self.soft_start = soft_start_stop[0] if not np.isnan(soft_start_stop[0]) else start
        self.soft_end = soft_start_stop[1] if not np.isnan(soft_start_stop[1]) else end

        # set soft start, soft end to be within start, end
        self.soft_start = min(max(self.soft_start, start), end)
        self.soft_end = max(min(self.soft_end, end), start)

        # do some checks
        assert self.window_size > 0, "Window size must be greater than 0"
        assert self.start < self.end, "Start must be less than end"
        assert self.strand in ['+', '-'], "Strand must be + or -"
        assert self.soft_start <= self.soft_end, "Soft start must be less than or equal to soft end"

        if direction not in ['head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref']:
            raise ValueError("Invalid direction, must be 'head-to-tail' or 'tail-to-head' or 'increasing-ref' or "
                             "'decreasing-ref'")

        # create windows
        self.windows = []

        if ((self.strand == '+' and self.direction == 'tail-to-head') or
                (self.strand == '-' and self.direction == 'head-to-tail') or
                self.direction == 'increasing-ref'):
            for k in range(self.start, self.end, self.window_size):
                self.windows.append((k, k + self.window_size))
        else:
            for k in range(self.end - 1, self.start - 1, -self.window_size):
                self.windows.append((k - self.window_size + 1, k + 1))

        # - loop through windows.
        # - if any window has only a partial overlap with the soft coordinates, then adjust it.
        # - if a window has no overlap with the bed entry, then it gets assigned a zero size automatically.
        self.windows = [(min(self.soft_end, max(self.soft_start, w[0])),
                         max(self.soft_start, min(self.soft_end, w[1]))) for w in self.windows]

    def __iter__(self):
        """ Return iterator """
        self.current_position = 0
        return self

    def __next__(self):
        """ Return next position """
        while self.current_position < len(self.windows):
            self.current_position += 1
            return self.windows[self.current_position - 1]
        raise StopIteration


class GenomeCoverage:
    def __init__(self, contigs: dict, stranded: bool = False) -> None:
        """
        Initialize the GenomeCoverage with contigs and stranded option.
        Written with ChatGPT.

        Args:
            contigs (dict): Dictionary with contig names as keys and their lengths as values.
            stranded (bool): Whether the coverage is stranded or not.

        Returns:
            None
        """
        self.contigs = contigs  # Dictionary with contig names and their lengths
        self.stranded = stranded  # Whether the coverage requested is stranded or not

        # Structure to store coverage data
        self.coverage = {}
        for contig, length in contigs.items():
            # check that contig is a string and length is an integer
            if not (isinstance(contig, str) and isinstance(length, int) and length > 0):
                raise ValueError("Contig names must be strings and lengths must be positive integers")
            if self.stranded:
                # If stranded, we keep two sets for each contig, one for positive and one for negative strand
                self.coverage[contig + "_+"] = list()
                self.coverage[contig + "_-"] = list()
            else:
                self.coverage[contig] = list()

        # Calculate total genome size
        self.total_bases = sum(self.contigs.values())

        # Double the genome size if stranded
        if self.stranded:
            self.total_bases *= 2

        # Initialize covered bases to 0
        self.covered_bases = 0

    def update_coverage_with_bed_line(self, bed_file_line: str) -> None:
        """  Update the coverage based on the given bed line.

        Args:
            bed_file_line (str): A line from the bed file.

        Returns:
            None
        """

        # Skip header lines
        if bed_file_line.startswith(("#", "track", "browser")):
            return

        columns = bed_file_line.strip().split("\t")
        if len(columns) < 3:  # BED line should have at least 3 columns (chrom, start, end)
            raise ValueError("Invalid BED line: {}".format(bed_file_line))

        chrom, start, end = columns[0], int(columns[1]), int(columns[2])
        strand = columns[5] if self.stranded else "."  # Assuming standard 6 column BED format

        # Check if chromosome is valid
        if chrom not in self.contigs:
            raise ValueError("Invalid chromosome: {}".format(chrom))

        # Check if start and end are valid
        if start < 0 or end > self.contigs[chrom] or start > end:
            raise ValueError("Invalid start or end position: {}".format(bed_file_line))

        # Record coverage
        if start < end:
            if self.stranded:
                if strand in ['+', '-']:
                    self.coverage[f"{chrom}_{strand}"].append([start, end])
            else:
                self.coverage[chrom].append([start, end])

    def is_genome_covered_within_tol(self, tolerance: int) -> bool:
        """ Check if genome is covered within the given tolerance.

        Args:
            tolerance (int): The tolerance value. So we check if genome size - covered region size <= tolerance.

        Returns:
            bool: True if genome is covered within the given tolerance, False otherwise.

        """
        # Handle a trivial case first
        if tolerance >= self.total_bases:
            return True

        # Calculate covered bases
        self.covered_bases = sum(get_length_of_merged_interval(covered) for contig, covered in self.coverage.items())

        if self.covered_bases > self.total_bases:
            raise ValueError("Something has gone wrong - covered base pairs cannot be greater than total base pairs")
        elif self.total_bases - self.covered_bases <= tolerance:
            return True
        else:
            return False

    def is_overlapping_intervals(self) -> bool:
        """ Check if there are overlapping intervals in the coverage.

        Returns:
            bool: True if there are overlapping intervals, False otherwise.

        """
        for contig, covered in self.coverage.items():
            covered_sorted = sorted(covered, key=lambda x: x[0])
            for i in range(len(covered_sorted) - 1):
                if covered_sorted[i][1] > covered_sorted[i + 1][0]:
                    return True
        return False


def intersect_two_genomic_intervals(start_1: int, end_1: int, strand_1: str,
                                    start_2: int, end_2: int, strand_2: str,
                                    norm_by_int_1_length: bool = False) -> Tuple[str, int, int]:
    """Intersect two genomic intervals and return the (signed) coordinates relative to the first interval

    Args:
        start_1 (int): start of the first interval
        end_1 (int): end of the first interval, must be >= start_1
        strand_1 (str): strand of the first interval, can be '+', '-' or '.'
        start_2 (int): start of the second interval
        end_2 (int): end of the second interval, must be >= start_2
        strand_2 (str): strand of the second interval, can be '+', '-' or '.'
        norm_by_int_1_length (bool): (default False) if True, normalize the coordinates by the length of first interval

    Returns:
        tuple: (type, coordinate_1, coordinate_2) where
        - type = 'same_strand' or 'opposite_strand' or 'unknown' or 'no_intersection'
        - coordinate_1 = start intersection of coordinate relative to the "arrow tail" of the first interval
        - coordinate_2 = end intersection of coordinate relative to the "arrow tail" of the first interval
        - By arrow tail, we mean start_1 if interval 1 is on the + strand or if strand_1 or strand_2 are '.',
          and end if interval 1 is on the - strand
        - coordinate_1 always <= coordinate_2.
        - the two coordinates can be positive or negative, depending on the orientation of the first interval
        - if no intersection, return ('no_intersection', 0, 0)

    Examples:
        >>> intersect_two_genomic_intervals(10, 20, '+', 15, 25, '+')
        ('same_strand', 5, 10)
        >>> intersect_two_genomic_intervals(10, 20, '+', 15, 25, '-')
        ('opposite_strand', 5, 10)
        >>> intersect_two_genomic_intervals(10, 20, '+', 15, 25, '.')
        ('unknown', 5, 10)
        >>> intersect_two_genomic_intervals(10, 20, '-', 15, 25, '+')
        ('opposite_strand', -5, 0)
        >>> intersect_two_genomic_intervals(10, 20, '-', 15, 25, '-')
        ('same_strand', -5, 0)
        >>> intersect_two_genomic_intervals(10, 20, '-', 15, 25, '.')
        ('unknown', 5, 10)
        >>> intersect_two_genomic_intervals(10, 20, '-', 15, 28, '.', norm_by_int_1_length=True)
        ('unknown', 0.5, 1.0)
        >>> intersect_two_genomic_intervals(10, 20, '-', 115, 125, '.')
        ('no_intersection', 0, 0)
    """

    combined_strand = strand_1 + strand_2

    # check if intervals are valid
    if start_1 >= end_1 or start_2 >= end_2:
        raise ValueError('start must be < end for both intervals')

    # set intersect type
    if combined_strand == '++' or combined_strand == '--':
        intersect_type = 'same_strand'
    elif combined_strand == '+-' or combined_strand == '-+':
        intersect_type = 'opposite_strand'
    elif combined_strand == '..' or combined_strand == '.+' or combined_strand == '.-' \
            or combined_strand == '+.' or combined_strand == '-.':
        intersect_type = 'unknown'
    else:
        raise ValueError('Invalid strand')

    # set arrow tail
    if (strand_1 == '.' or strand_2 == '.') or strand_1 == '+':
        arrow_tail = start_1
    elif strand_1 == '-':
        arrow_tail = end_1
    else:
        raise ValueError('Invalid strand')

    # calculate coordinates
    coordinate_1 = max(start_1, start_2) - arrow_tail
    coordinate_2 = min(end_1, end_2) - arrow_tail

    # if norm_by_int_1_length, normalize coordinates by length of first interval
    if norm_by_int_1_length:
        coordinate_1 /= (end_1 - start_1)
        coordinate_2 /= (end_1 - start_1)

    if coordinate_1 >= coordinate_2:
        return 'no_intersection', 0, 0
    else:
        return intersect_type, coordinate_1, coordinate_2


def get_length_of_merged_interval(interval_list: list[list[int, int]]) -> int:
    """ Merge the intervals in a list and return the total length of the merged interval.

    Args:
        interval_list: A list of lists, each containing two integers representing the start and end of an interval.

    Returns:
        int: The total length of the merged interval.

    Examples:
        >>> get_length_of_merged_interval([[1, 3], [2, 6], [8, 10], [15, 18]])
        10
        >>> get_length_of_merged_interval([[1, 4], [4, 5]])
        4
    """
    if len(interval_list) == 0:
        return 0

    # copy interval_list to avoid modifying the original list
    interval_list_copy = copy.deepcopy(interval_list)

    # Sort the intervals by their start position
    interval_list_copy.sort(key=lambda x: x[0])

    # Merge the intervals
    merged_intervals = [interval_list_copy[0]]
    for interval in interval_list_copy[1:]:
        if interval[0] <= merged_intervals[-1][1]:
            # Overlap with the previous interval, merge them
            merged_intervals[-1][1] = max(merged_intervals[-1][1], interval[1])
        else:
            # No overlap, add the interval to the list
            merged_intervals.append(interval)

    # Calculate the total length of the merged interval
    total_length = 0
    for interval in merged_intervals:
        total_length += interval[1] - interval[0]

    return total_length
