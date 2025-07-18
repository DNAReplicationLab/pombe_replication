import sys
from validate_bed_format import FilteredBedFileIterator


def add_fork_boundaries_to_rDNA_combined_bed_line(bed_line: str) -> str:
    """ Create a detectIndex column for a line from rDNA combined bed that contains fork boundaries
    Args:
        bed_line: Tab-separated line from a bed file.
                  Using 0-indexing, the seventh column is the step, the tenth and eleventh columns are the left
                  and right boundaries of the interval over which the step was identified, and the last column
                  is the detectIndex of format "readID_contig_start_end_orientation".

    Returns:
        The input line with the detectIndex column updated to include the boundaries of the interval over which the
         step was identified.

    """

    columns = bed_line.strip().split('\t')
    step = float(columns[7])
    boundary_left_of_step = columns[10]
    boundary_right_of_step = columns[11]

    if step < 0:
        detectIndex = f"_L_{boundary_left_of_step}_{boundary_right_of_step}"
    else:
        detectIndex = f"_R_{boundary_left_of_step}_{boundary_right_of_step}"

    columns[-1] = columns[-1] + detectIndex

    return '\t'.join(columns)


if __name__ == '__main__':

    # check that inputs have been piped in
    if sys.stdin.isatty():
        print('Usage: cat <file> | python rDNA_process_calc_steps_into_pause_file.py > <output_file>')
        sys.exit(1)
    
    # Make and print the header
    header_list = ["contig", "pauseSite", "pauseSitePlus1", "read_id", "abs_step_normalized_by_sd", "strand",
                             "alignLen", "step", "dens_at_left_of_step", "dens_at_right_of_step",
                             "boundary_left_of_step", "boundary_right_of_step", "detectIndex"]
    header_line = "\t".join(header_list)
    print(header_line)

    # Read input from stdin and process each line
    for line in FilteredBedFileIterator(sys.stdin, allow_float_score=True):

        # check that the length of the line matches the length of the header list
        if len(line.strip().split('\t')) != len(header_list):
            raise ValueError(f"ERROR: Line does not match header length. Line: {line}")

        # Process the line and print the result
        print(add_fork_boundaries_to_rDNA_combined_bed_line(line))
