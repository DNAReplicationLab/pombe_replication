import sys
from DNAscentTools.modBAM_tools_additional import process_fork_index

if __name__ == "__main__":

    desc = """    Filter pause-bed-file lines to remove intervals that lie outside forks/alignments
    
    Logic: "Illegal" intervals are those that lie outside the fork/alignment, so they are removed.
           Even those that lie only partially outside the fork/alignment are removed. 

    Input: - Pipe in a bed file that was obtained by converting a pause file to a bed file.
           - Such a file can be obtained by running the script convert_pause_file_to_bed.py and will have the following
             format:
             * The first six columns are the standard BED6 columns, and the rest are the pause file columns.
               The standard BED6 columns are: contig, start, end, readID, score, orientation.
             * In the remaining columns, exactly one column is the detectIndex corresponding to the fork
               and will have the format "readID_contig_startAlign_endAlign_orn_dirn_startFork_endFork".
               The labels in the string are self-explanatory, except orn=fwd/rev and dirn=L/R.
               startFork always < endFork irrespective of whether it's an L or a R fork.
           - We check if the interval lies within the fork/alignment using the coordinates in the detectIndex.
           - If there are multiple columns that can be interpreted as detectIndices,
             the interval is rejected if it lies outside any of the forks/alignments.
           - We do not expect multiple detectIndices to be present, but if they are, they are used as described above.
           - If there are no detectIndices on a line, that interval is always accepted.

    Output: Same lines as input, except illegal intervals are removed.
    
    Example of logic using a sample bed line (tabs have been replaced by spaces for readability):
    chr1 10000 12000 c773aee5-3f41-4ef4-a0f7-c5f6fbc5548f 0 + c773aee5-3f41-4ef4-a0f7-c5f6fbc5548f_chr1_0_40000_fwd_L_11000_22000
    From the line above, it is clear that the given read is aligned from 0 to 40000 on the forward strand of chr1 and
    has a left fork in the interval 11000-22000. So, if the user invokes filtration for forks, the interval 10000-12000
    lies outside the fork 11000-22000, so it will be rejected. But if the user invokes filtration for alignments, the
    interval 10000-12000 lies within the alignment 0-40000, so it will be accepted.

    Sample usage:
        < pauseFile.bed python filter_pause_file_bed_intervals_outside.py --forks > pauseFile_filtered.bed 
        < pauseFile.bed python filter_pause_file_bed_intervals_outside.py --alignments > pauseFile_filtered.bed
    """

    # check that the input is coming from a pipe
    if sys.stdin.isatty():
        print(desc)
        sys.exit(1)

    is_filter_forks = False
    is_filter_alignments = False
    # read the first argument to decide whether to filter forks or alignments
    if len(sys.argv) != 2:
        print("Please provide exactly one argument: --forks or --alignments")
        sys.exit(1)
    elif sys.argv[1] == "--forks":
        is_filter_forks = True
        is_filter_alignments = False
    elif sys.argv[1] == "--alignments":
        is_filter_forks = False
        is_filter_alignments = True
    else:
        print("Please provide exactly one argument: --forks or --alignments")
        sys.exit(1)

    for line in sys.stdin:
        fields = line.strip().split(sep="\t")

        # print comment lines and bed header lines as is
        if fields[0].startswith(("#", "track", "browser")):
            print(line)
            continue

        # get the standard bed fields
        bed_contig = fields[0]
        bed_start = int(fields[1])
        bed_end = int(fields[2])

        # set the flag to reject the entry
        reject_entry = False

        # iterate over the remaining fields to find the detectIndex
        for field in fields[3:]:
            try:
                read_id, contig, start, end, strand, fork_type, fork_start, fork_end = process_fork_index(field)
                # NOTE: even if the bed file contains strand information, we do not check if the strand matches
                #       between the bed entry and the fork index because often we set strands in the bed file using the
                #       fork direction or whether the fork is leading or lagging to perform some calculations
                #       and this is not expected to match the strand in the fork index.
                if ((is_filter_forks and (bed_start < fork_start or bed_end > fork_end)) or
                   (is_filter_alignments and (bed_start < start or bed_end > end)) or
                   (bed_contig != contig) or (bed_start > bed_end)):
                    reject_entry = True
                    break
            except ValueError:
                continue
            except IndexError:
                continue

        # print the line if it's not rejected
        if not reject_entry:
            print(line, end="")
