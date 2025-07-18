import sys
from perform_binary_bedgraph_operation import process_bedgraphs


if __name__ == "__main__":

    # check inputs
    if len(sys.argv) < 3:
        print("Usage: python correlation_coefficient_bedgraph_values.py <input_file_1> <input_file_2>")
        print("Goal: Given 2 bedgraph files with identical coordinates, compute the correlation coefficient of values.")
        print("      Mathematically speaking, calculate https://en.wikipedia.org/wiki/Pearson_correlation_coefficient")
        print("      where the two vectors are the values in the two input files.")
        print("Error: Please provide exactly two input files.")
        sys.exit(1)

    # perform calculation
    # NOTE: we do the calculation in a slightly clunky way as we already have functions for the basic
    #       bedgraph operations.
    df_multiply_a_b = process_bedgraphs(sys.argv[1], sys.argv[2], "multiply")
    num_values = df_multiply_a_b.shape[0]
    df_sum_a_b = df_multiply_a_b["result"].sum()

    df_twice_a = process_bedgraphs(sys.argv[1], sys.argv[1], "add")["result"]
    df_twice_b = process_bedgraphs(sys.argv[2], sys.argv[2], "add")["result"]
    num_values_a = df_twice_a.shape[0]
    num_values_b = df_twice_b.shape[0]
    df_sum_a = df_twice_a.sum() / 2
    df_sum_b = df_twice_b.sum() / 2

    df_sum_a_sq = process_bedgraphs(sys.argv[1], sys.argv[1], "multiply")["result"].sum()
    df_sum_b_sq = process_bedgraphs(sys.argv[2], sys.argv[2], "multiply")["result"].sum()

    if num_values > 0:
        numerator = df_sum_a_b - (df_sum_a * df_sum_b) / num_values
        denominator = ((df_sum_a_sq - (df_sum_a ** 2) / num_values) *
                       (df_sum_b_sq - (df_sum_b ** 2) / num_values)) ** 0.5
        if denominator == 0 or not (num_values == num_values_a == num_values_b):
            print("NaN")
        else:
            print(f"{numerator / denominator}")
