import pandas as pd
import sys

if __name__ == "__main__":

    desc = """ Goal: Given bedgraph with a signal, calculate its spatial correlation.
    
    Usage: < bedgraph python calculate_bedgraph_spatial_autocorrelation.py \
      max_distance > statistics.txt
      
    Input: * Bedgraph with some signal (e.g. sensitivity) in the 4th column. 
           * Bedgraph must have all intervals along the genome and each interval must have the same length.
           * Max distance is the maximum distance to consider for calculation, must be a multiple of the
             bedgraph interval length, or it will be rounded down to the nearest multiple. 
             
    Output: * Two tab-separated columns: lag_window_size, correlation
            * Comments start with # and column names are present.
            * We are calculating < (s(x) - <s>) * (s(x+lag) - <s>) > / (<s>)^2 where s is the
              signal, lag is the lag window size, and <> is the average over all x
              (after grouping by contig, to ensure we don't mix contigs).
      CAUTION: * You should set max_distance as a value that is well into the tail otherwise
                 you will get truncated results.
               * We recommended something in the range of a few kb to tens of kb, basically a few times the mean
                 distance between neighbouring pauses.
               * Do not set this to a large value, then the calculation will take a very long time and a lot of
                 this will be wasted time calculating probabilities that pauses are far away from each other.
    """

    # check if input is piped in, otherwise print description and exit
    if sys.stdin.isatty():
        print(desc)
        sys.exit(1)

    # read in bedgraph to a pandas dataframe
    df = pd.read_csv(sys.stdin, sep=" ", header=None, names=["chr", "start", "end", "signal"], comment="#")

    # get interval length
    df["interval_length"] = df["end"] - df["start"]
    interval_length = df["interval_length"].median()

    # check that intervals are of positive length
    if df["interval_length"].min() <= 0:
        print("ERROR: all intervals must have positive length, exiting.")
        sys.exit(1)

    # check that interval lengths don't vary by much
    interval_length_mean = df["interval_length"].mean()
    interval_length_std = df["interval_length"].std()
    if interval_length_std > interval_length_mean / 20:
        print("ERROR: we need uniform interval lengths for this calculation.")
        print("Interval lengths sd is more than 5% of the mean, exiting.")
        sys.exit(1)

    # calculate mean and sd of signal over all intervals
    mean_signal = df["signal"].mean()
    sd_signal = df["signal"].std()

    # subtract mean signal from signal
    df["sigMinusMean"] = df["signal"] - mean_signal

    # read max distance from command line, if not given, exit
    try:
        max_distance = int(sys.argv[1])
        if max_distance < 1:
            raise ValueError
    except IndexError:
        print("No max distance given, exiting.")
        sys.exit(1)
    except ValueError:
        print("Max distance must be an integer, exiting.")
        sys.exit(1)

    # round max distance down to nearest multiple of interval length
    max_distance = max_distance - (max_distance % interval_length)

    # get number of bins
    num_bins = int(max_distance / interval_length)

    # group df by chr
    df_grouped = df.groupby("chr")

    # create two lists to store lag window size and correlation
    lag_window_size = []
    correlation = []

    # loop through bins
    for i in range(1, num_bins):

        lag_window_size.append((i - 1) * interval_length)

        correlation_grouped = 0
        count = 0

        # loop through chromosomes
        for chr_name, chr_df in df_grouped:

            # get the sum of the product of all entries in the window and count the number of windows
            correlation_grouped += chr_df["sigMinusMean"].rolling(i).apply(lambda s: s.iloc[0] * s.iloc[i - 1]).sum()
            count += chr_df.shape[0] - i + 1

        correlation.append(correlation_grouped/count)

    # normalize pdf values to 1
    correlation_normalized = [x / (mean_signal ** 2) for x in correlation]

    # create dataframe from lists
    df_out = pd.DataFrame({"lag_window_size": lag_window_size, "correlation": correlation_normalized})

    # print mean, sd values to stdout as a comment
    print("# mean signal: " + str(mean_signal))
    print("# sd signal: " + str(sd_signal))
    print("# correlation normalized to mean signal squared")
    
    # print dataframe to stdout
    df_out.to_csv(sys.stdout, sep="\t", index=False)
