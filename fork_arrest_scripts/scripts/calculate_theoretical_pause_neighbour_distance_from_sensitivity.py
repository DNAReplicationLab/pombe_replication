import pandas as pd
import sys
import numpy as np

if __name__ == "__main__":

    desc = """ Goal: Given calculated pause sensitivity bedgraph, calculate theoretical distance statistics between
    neighbouring pauses.
    
    Usage: < pause_sensitivity.bedgraph python calculate_theoretical_pause_neighbour_distance_from_sensitivity.py \
      max_distance [bin_size] > distance_statistics.txt
      
    Input: * Bedgraph with calculated pause sensitivity scores and max distance to consider between neighbouring pauses. 
           * Bedgraph must have all intervals along the genome and each interval must have the same length and all
             entries must be between 0 and 1.
           * Bedgraph values must sum up to the total number of pauses (this will be true if this was generated using
             our sensitivity-calculation scripts).
           * Max distance is the maximum distance to consider between neighbouring pauses, must be a multiple of the
             bedgraph interval length, or it will be rounded down to the nearest multiple.
           * Bin size (an optional parameter) is the size of the bins to use in the calculation, defaults to the
             interval length of the bedgraph. If bin size is a multiple of the bedgraph interval length,
             we will sum up the bedgraph intervals to get appropriate intervals.
             Otherwise, the script will throw an error. 
           
    Mathematical logic: 
           1. Given a point of interest P, a grid of some size, and a distance k,
             the probability that the nearest pause to P is at distance k =
             (probability of no pause at P + 1) * (probability of no pause at P - 1) * \
             (probability of no pause at P + 2) * (probability of no pause at P - 2) * ... \
             (probability of no pause at P + k - 1) * (probability of no pause at P - (k-1)) * ... \
             ((probability of pause at P + k) + (probability of pause at P - k))
           2. Calculate the product in the step above for all points P and sum/average.
           3. Repeat step 2 for all distances k from grid size to max distance.
           4. Normalize the sum/average from step 3 so that the sum of all probabilities is equal to the
              detected number of pauses. This is so that we can compare the theoretical distance statistics
              directly with the detected distance statistics.
           NOTE: the probability of pause at any point P is equal to the sensitivity at P from the bedgraph.
             
    Output: * Eight, tab-separated columns: bin_lower_limit, bin_higher_limit, theoretical_count, theoretical_sd,
                                            uniform_theoretical_count, uniform_theoretical_sd, pdf, uniform_pdf
            * Comments start with # and column names are present.
            * The "uniform" columns are calculated by assuming that the pauses are uniformly distributed
              along the genome.
            * The pdf columns are a probability distribution function and add up to 1.
            * The count columns are the pdf multiplied by the total number of pauses.
            * The sd columns are the standard deviation of the pdf, calculated using standard counting statistics.
      CAUTION: * You should set max_distance as a value that is well into the tail of the histogram otherwise
                 you will get truncated results.
               * We recommended something in the range of a few kb to tens of kb, basically a few times the mean
                 distance between neighbouring pauses.
               * Do not set this to a large value, then the calculation will take a very long time and a lot of
                 this will be wasted time calculating probabilities that pauses are far away from each other.
               * If you expect to see a length scale of order the contig length in the distance statistics,
                 then something is going wrong.
      CAUTION 2: If the count is much larger than the sd in a bin, i.e. count - sd is well away from zero, then one
                 can interpret the distribution in that bin as gaussian. So, one can ask questions like 3 sigma away
                 etc. with the usual understanding of what that means. 
                 However, if the count - sd is small enough to be close to zero, then the distribution is not gaussian
                 and 3 sd may not mean 99.7% confidence etc. One will run into such situations in the tail.
                 This is not a problem with the program, but something the user should be aware of. 
    """

    # check if input is piped in, otherwise print description and exit
    if sys.stdin.isatty():
        print(desc)
        sys.exit(1)

    # read in bedgraph to a pandas dataframe
    df = pd.read_csv(sys.stdin, sep=" ", header=None, names=["chr", "start", "end", "sensitivity"], comment="#")

    # sort df by chr and start
    df = df.sort_values(by=["chr", "start"])

    # get interval length
    df["interval_length"] = df["end"] - df["start"]
    interval_length = df["interval_length"].median()

    # check that interval lengths don't vary by much
    interval_length_mean = df["interval_length"].mean()
    interval_length_std = df["interval_length"].std()
    if interval_length_std > interval_length_mean / 20:
        print("ERROR: we need uniform interval lengths for this calculation.")
        print("Interval lengths sd is more than 5% of the mean, exiting.")
        sys.exit(1)

    # calculate total number of pauses, total length of genome, and half of the mean separation between pauses
    N_total_pauses = df["sensitivity"].sum()
    L_genome = df["interval_length"].sum()
    mean_separation_pauses = L_genome / N_total_pauses

    # drop the interval length, start, and end columns
    df = df.drop(columns=["interval_length", "start", "end"])

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

    # try reading bin size from command line, if not given, set it to interval length
    try:
        bin_size = int(sys.argv[2])

        assert bin_size >= interval_length
        assert bin_size % interval_length == 0
        assert bin_size < mean_separation_pauses / 7

        # if bin size is greater than interval length, sum up the bedgraph intervals to get appropriate intervals
        if bin_size > interval_length:
            group = df.groupby('chr').cumcount() // int(bin_size / interval_length)
            df = df.groupby(['chr', group]).agg({'sensitivity': 'sum'}).droplevel(1).reset_index()
            interval_length = bin_size

    except IndexError:
        bin_size = interval_length
    except AssertionError:
        print("Problems with bin size, must be a multiple of interval length and >= interval length, exiting.")
        sys.exit(1)
    except BaseException:
        print("Something went wrong with bin size, exiting.")
        sys.exit(1)

    # round max distance down to nearest multiple of bin size
    max_distance = max_distance - (max_distance % bin_size)

    # get number of bins
    num_bins = int(max_distance / bin_size)

    # ensure there are at least a few bins, like 10
    if num_bins < 10:
        print("Max distance is too small, we need at least 10 bins, exiting.")
        sys.exit(1)

    # create three lists to store bin lower limit, bin higher limit, and pdf value
    bin_lower_limit = []
    bin_higher_limit = []
    pdf_value = []

    # NOTE: in these calculations, we are assuming that the mean separation between pauses is much smaller
    # than the length of any contig. If this is not the case, the theoretical expectation from uniform
    # distribution could be off. So, we check this below.
    min_contig_length = df.groupby("chr")["sensitivity"].count().min() * interval_length
    assert(min_contig_length > 10 * mean_separation_pauses)

    # ensure all sensitivity values are between 0 and 1
    if any(df["sensitivity"] < 0) or any(df["sensitivity"] > 1):
        print("Sensitivity values must be between 0 and 1, exiting.")
        sys.exit(1)

    # create a column 1 - sensitivity
    df["oneMinusSensitivity"] = 1 - df["sensitivity"]

    # loop through bins
    for i in range(1, num_bins + 1):
        bin_lower_limit.append((i - 1) * bin_size)
        bin_higher_limit.append(i * bin_size)

        # create a window of size 2 * i + 1 around each entry in chr_df and calculate the product
        # (1-p_1) * (1-p_(-1)) * (1-p_2) * (1-p_(-2)) * ... * (1-p_(i-1)) * (1-p_(-i+1)) * (p_i + p_(-i))
        # where position 0 is the current entry, and sum all products
        # NOTE: since we've already calculated 1 - sensitivity, the product looks somewhat different in the
        # code below, but it is the product written above
        def calculate_special_product(x):
            # if any entry in x is nan or length of x is 0, return 0
            if len(x) == 0 or any(x.isna()):
                return 0
            else:
                N = len(x) // 2
                indices_to_consider = list(range(1, N)) + list(range(N + 1, 2 * N))
                return (2 - x.iloc[2*N] - x.iloc[0]) * x.iloc[indices_to_consider].prod()
                # NOTE: when indices_to_consider is empty, x.iloc[indices_to_consider].prod() returns 1

        # get the sum of the product of all entries in the window
        pdf_value.append(df.groupby("chr")["oneMinusSensitivity"].rolling(2 * i + 1).
                         apply(calculate_special_product).mean())

    # if any pdf_value is negative or zero, exit
    if any([x <= 0 for x in pdf_value]):
        print("PDF value is negative or zero, exiting.")
        sys.exit(1)

    # using pdf, calculate theoretical value and sd
    theoretical_value = [x * N_total_pauses for x in pdf_value]
    theoretical_sd = [np.sqrt(x * (1 - x) * N_total_pauses) for x in pdf_value]

    # calculate theoretical value from random distribution
    uniform_calc = [np.exp(-2 * k[0] / mean_separation_pauses) - np.exp(-2 * k[1] / mean_separation_pauses)
                    for k in zip(bin_lower_limit, bin_higher_limit)]

    # calculate theoretical value and sd from random distribution
    uniform_theoretical_value = [x * N_total_pauses for x in uniform_calc]
    uniform_theoretical_sd = [np.sqrt(x * (1 - x) * N_total_pauses) for x in uniform_calc]

    # create dataframe from lists
    df_out = pd.DataFrame({"bin_lower_limit": bin_lower_limit, "bin_higher_limit": bin_higher_limit,
                           "theoretical_count": theoretical_value, "theoretical_sd": theoretical_sd,
                           "uniform_theoretical_count": uniform_theoretical_value,
                           "uniform_theoretical_sd": uniform_theoretical_sd,
                           "pdf": pdf_value,
                           "uniform_pdf": uniform_calc})

    # print dataframe to stdout
    df_out.to_csv(sys.stdout, sep="\t", index=False)
