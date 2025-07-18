import argparse
import numpy as np
import datetime

if __name__ == "__main__":

    desc = """    Simulate many step functions, fit a step separately to each, calculate error and other statistics.

    Input: No input file required. Just supply the necessary parameters.

    Output: Two files are output: a list of the errors in the step detection, and the aggregate step profile
            after aligning the simulated curves according to the detected steps.

    Sample usage:
        python calculate_theoretical_error_in_step_detection.py --num_simulations 1000 --lo 0.1 --hi 0.9 \
          --num_low_points 100 --num_high_points 100 --output_file_prefix test --trials_per_point 100
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--num_simulations', type=int, required=True,
                        help="Number of simulations to run.")
    parser.add_argument('--lo', type=float, required=True,
                        help="Low value of the step.")
    parser.add_argument('--hi', type=float, required=True,
                        help="High value of the step.")
    parser.add_argument('--num_low_points', type=int, required=True,
                        help="Number of data points in the low value of the step.")
    parser.add_argument('--num_high_points', type=int, required=True,
                        help="Number of data points in the high value of the step.")
    parser.add_argument('--trials_per_point', type=int, required=True,
                        help="Number of trials per point to use in the Binomial distribution.")
    parser.add_argument('--output_file_prefix', type=str, required=True,
                        help="Prefix for the output files. If the prefix is 'test', the output files will be "
                             "'test_histogram.png' and 'test_aggregate_step_profile.png'. If you want files to go "
                             "to a different directory, you can specify the full path here. e.g. "
                             "'/home/user/test' will output files to '/home/user/test_histogram.png' and "
                             "'/home/user/test_aggregate_step_profile.png'.")

    args = parser.parse_args()

    # extract the arguments
    n_simulations = args.num_simulations
    lo = args.lo
    hi = args.hi
    n_low = args.num_low_points
    n_high = args.num_high_points
    output_file_prefix = args.output_file_prefix
    trials_per_point = args.trials_per_point

    # string with parameters to be saved in the output file
    parameter_str = ("# { \"num_simulations\": %d, \"lo\": %.3f, \"hi\": %.3f, "
                     "\"num_low_points\": %d, \"num_high_points\": %d, "
                     "\"trials_per_point\": %d, \"date\": %s }\n" % (n_simulations, lo, hi, n_low, n_high,
                                                                     trials_per_point, str(datetime.datetime.now())))

    # create empty list to store the error in the step position
    step_error = []

    # initialize object to store all simulated curves
    simulated_data_list = []

    # perform simulations
    for i in range(n_simulations):

        # simulate a step function
        # Draw n_high points from a binomial distribution with probability hi and trials_per_point trials,
        # and n_low points from a binomial distribution with probability lo and trials_per_point trials.
        y = (list(np.random.binomial(n=trials_per_point, p=hi, size=n_high) / trials_per_point) +
             list(np.random.binomial(n=trials_per_point, p=lo, size=n_low) / trials_per_point))

        # initialize the known step-point and the minimum goodness of fit
        step_ground_truth = n_high
        gof_min = 1e10  # a large number

        # initialize the best-fit step location
        step_best_fit = 0

        # fit a step to the simulated data and find location of the step
        for location in range(max(0, step_ground_truth - 20), min(n_low + n_high, step_ground_truth + 20 + 1)):

            y_hat = [hi] * location + [lo] * (n_low + n_high - location)
            residuals = [y[i] - y_hat[i] for i in range(n_low + n_high)]
            gof = sum([k ** 2 for k in residuals])

            if gof < gof_min:
                gof_min = gof
                step_best_fit = location

        # store error in the step position
        step_error.append(step_best_fit - step_ground_truth)

        # store the data
        simulated_data_list.append(y)

    # average the simulated curves after aligning by the detected steps
    mean_step_profile = [0] * (n_low + n_high)
    counts = [0] * (n_low + n_high)

    for i in range(n_simulations):
        for j in filter(lambda x: 0 <= x - step_error[i] < n_low + n_high, range(n_low + n_high)):
            mean_step_profile[j] += simulated_data_list[i][j - step_error[i]]
            counts[j] += 1

    mean_step_profile = [mean_step_profile[i] / counts[i] for i in range(n_low + n_high) if counts[i] > 0]

    # calculate sd of the step profile after aligning by the detected steps
    sd_step_profile = [0] * (n_low + n_high)

    for i in range(n_simulations):
        for j in filter(lambda x: 0 <= x - step_error[i] < n_low + n_high, range(n_low + n_high)):
            sd_step_profile[j] += (simulated_data_list[i][j - step_error[i]] - mean_step_profile[j]) ** 2

    sd_step_profile = [np.sqrt(sd_step_profile[i] / counts[i]) for i in range(n_low + n_high) if counts[i] > 0]

    # save the error in the step position to a file
    with open(output_file_prefix + "list_step_error.txt", "w") as f:

        # sort the list of step errors
        step_error.sort()

        # write the list of step errors to a file
        f.write(parameter_str)
        f.write("step_error\n")
        for item in step_error:
            f.write("%d\n" % item)

    # save the aggregate step profile to a file
    with open(output_file_prefix + "aggregate_step_profile.txt", "w") as f:
        f.write(parameter_str)
        f.write("mean\tsd\tcount\n")
        for item in zip(mean_step_profile, sd_step_profile, counts):
            f.write("%.3f\t%.3f\t%d\n" % item)
