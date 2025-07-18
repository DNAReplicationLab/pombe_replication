import argparse
import sys
import os
from DNAscentTools.modBAM_tools_additional import ModBamRecordProcessor
from itertools import count

if __name__ == "__main__":

    desc = """    Output autocorrelation (over a window size) per read id from mod bam file. 

    Input: Pipe in modBAM file contents in plain text including or excluding headers.
    
    Logic: For each read in mod BAM, threshold each base and chop into (default non-overlapping) chunks of so many
           thymidines per chunk. Then, smooth the data according to the given window size (non-overlapping) and
           calculate autocorrelation for each chunk. 

    Output: Four columns of tab-separated data are output to stdout and and two text files are created.
        Stdout:
        First line is 'detectIndex chunk_index window_index auto_correlation_value'
        Every other line is a detectIndex and two numbers are given by the first line.
        detectIndex is a string in the format readID_contig_start_end_orn or just the readID (depending on whether
        a custom XA tag is present in the modBAM file).
        contig, start, end are positions on the reference genome and orn is fwd/rev.
        Text file 1:
        The first line is 'detectIndex chunk_index auto_correlation_length'
        So for each chunk, we calculate at which window coordinate does autocorrelation drop below 0.5
        and report it here (0.5 is the default value which can be changed through input options).
        Text file 2:
        The first line is 'detectIndex max_auto_correlation_length'
        After the calculation for file 1 above, we calculate the maximum value of auto_correlation_length
        per detectIndex and report it here.
        NOTE: Chunks with adjacent indices should not be interpreted as chunks that are adjacent in the read
        as we may drop some chunks due to different reasons but number them sequentially.
        NOTE: for best results, use a window size that is much smaller than the chunk size e.g. 300 and 10000.
        

    Sample usage:
        samtools view sample.mod.bam | python <programName.py> --thres 0.5 --window 300 --chunk 10000 --tag T
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--thres', type=float, required=False,
                        help='(default 0.5) threshold above (below) which T is regarded as modified (unmodified)',
                        default=0.5)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')
    parser.add_argument('--window', type=int, required=False,
                        help='(default 300) (non-overlapping) window size in thymidines to smooth the data',
                        default=300)
    parser.add_argument('--chunk', type=int, required=False,
                        help='(default 10000) number of thymidines per chunk',
                        default=10000)
    parser.add_argument('--chunkOverlap', type=int, required=False,
                        help='(default 0) number of thymidines to overlap between chunks',
                        default=0)
    parser.add_argument('--autoCorrelationThres', type=float, required=False,
                        help='(default 0.5) threshold used to measure autocorrelation length '
                             'i.e. at which window coordinate does autocorrelation drop below this value. '
                             'Must be between 0 and 1',
                        default=0.5)
    parser.add_argument('--textFilePrefix', type=str, required=True,
                        help='prefix for text files to be created e.g. if you want files named '
                             '/A/B/CD_summary_chunk.txt and /A/B/CD_summary_read.txt, '
                             'provide /A/B/CD here; the program will create any required directories')

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # check that window, chunk, chunkOverlap are positive and make sense
    if args.window <= 0 or args.window > args.chunk:
        raise ValueError("Bad window size. Must be positive and less than chunk size.")
    if args.chunk <= 0 or args.chunk <= args.chunkOverlap:
        raise ValueError("Bad chunk size. Must be positive and greater than chunk overlap.")

    # check that autoCorrelationThres is between 0 and 1
    if args.autoCorrelationThres < 0 or args.autoCorrelationThres > 1:
        raise ValueError("Bad autoCorrelationThres. Must be between 0 and 1.")

    # make any required directories
    # get parent directory of prefix
    parent_dir = os.path.dirname(args.textFilePrefix)

    # if parent directory doesn't exist, make it
    if not (parent_dir == "" or os.path.exists(parent_dir)):
        os.makedirs(parent_dir, exist_ok=True)

    # initialize variables
    modbam_record = ModBamRecordProcessor(args.thres, args.tag, True, True)

    # output data
    # first, initialize headers and files
    with open(f"{args.textFilePrefix}_summary_chunk.txt", "w") as fp_chunk:
        with open(f"{args.textFilePrefix}_summary_read.txt", "w") as fp_read:

            fp_chunk.write("detectIndex\tchunk_index\tauto_correlation_length\n")
            fp_read.write("detectIndex\tmax_auto_correlation_length\n")

            print("detectIndex\tchunk_index\twindow_index\tauto_correlation_value")
            for line in sys.stdin:
                if line.startswith("@"):
                    continue
                else:
                    window_index_max = 0
                    modbam_record.process_modbam_line(line.rstrip())
                    for chunk_index, k in zip(count(),
                                              modbam_record.calculate_autocorrelations(args.chunk, args.window,
                                                                                       args.chunkOverlap/args.chunk)):
                        is_threshold_achieved = False
                        for window_index, autocorr in zip(count(), k):

                            if autocorr <= args.autoCorrelationThres and not is_threshold_achieved:
                                fp_chunk.write(f"{modbam_record.read_id}\t{chunk_index}\t{window_index}\n")
                                is_threshold_achieved = True
                                if window_index > window_index_max:
                                    window_index_max = window_index

                            print(f"{modbam_record.read_id}\t{chunk_index}\t{window_index}\t{autocorr}")

                    if window_index_max > 0:
                        fp_read.write(f"{modbam_record.read_id}\t{window_index_max}\n")
