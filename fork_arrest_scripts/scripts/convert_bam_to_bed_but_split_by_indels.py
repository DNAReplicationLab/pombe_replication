import argparse
import sys
from DNAscentTools.modBAM_tools_additional import cigar_to_ref_to_query_tbl

if __name__ == "__main__":

    desc = """    Convert BAM to BED but reads are split if insertions/deletions beyond given threshold(s) are present.

    Input: Pipe in BAM or modBAM file contents in plain text including or excluding headers.

    Output: Six columns are output to stdout separated by tabs and with no column names:
        contig, start, end, read_id, score, strand. These are the standard BED columns, with score
        set to 1000 if primary alignment, and 0 if otherwise. Every read in the BAM file is output
        only once if no insertions or deletions are present. If insertions or deletions are present,
        the read is split into multiple lines, one for each segment between insertions or deletions.
        Unmapped reads are not output.

    Sample usage:
        samtools view sample.mod.bam | python <programName.py> --inThreshold 1000 --delThreshold 1000
            # split reads at insertions or deletions of 1000 bp or more
    """

    # get options
    large_num = 1_000_000_000_000_000
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inThreshold', type=int, required=False, default=large_num,
                        help='(default 1_000_000_000_000_000 i.e. unused) Threshold for insertions to split reads.')
    parser.add_argument('--delThreshold', type=int, required=False, default=large_num,
                        help='(default 1_000_000_000_000_000 i.e. unused) Threshold for deletions to split reads.')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # calculate and output data
    for line in sys.stdin:
        if line.startswith("@"):
            # this is a header line
            continue
        else:
            # process record and examine mapped reads
            read_id, flag, contig, start, quality, cigar, _ = line.rstrip().split("\t", maxsplit=6)
            start = int(start) - 1  # convert to 0-based start
            flag = int(flag)

            if not flag & 4 == 4:

                # set score and strand
                score = 0 if (flag & 256 == 256 or flag & 2048 == 2048) else 1000
                strand = '-' if flag & 16 == 16 else '+'

                # get table of mapping from reference to query, and get indel positions
                # NOTE: we ignore the first and last positions in the table as there are either softclips/hardclips
                #       or matches there. We do not expect reads to start or end with indels.
                tbl = cigar_to_ref_to_query_tbl(cigar, start)
                indel_positions = list(filter(lambda x: (tbl[x + 1][0] == tbl[x][0] and
                                                         tbl[x + 1][1] - tbl[x][1] >= args.inThreshold) or
                                                        (tbl[x + 1][1] == tbl[x][1] and
                                                         tbl[x + 1][0] - tbl[x][0] >= args.delThreshold),
                                              range(1, len(tbl) - 2)))

                # prepare start and end of bed intervals
                positions = [start] + [tbl[k][0] for k in indel_positions] + \
                            [tbl[k + 1][0] for k in indel_positions] + [tbl[-1][0]]
                positions.sort()

                # output bed intervals
                for k in zip(positions[::2], positions[1::2]):
                    print(f"{contig}\t{k[0]}\t{k[1]}\t"
                          f"{read_id}\t{score}\t{strand}")
