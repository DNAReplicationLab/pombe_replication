import sys
from perform_binary_bedgraph_operation import main


if __name__ == "__main__":

    # check inputs
    if len(sys.argv) < 3:
        print("Usage: python divide_bedgraph_values.py <input_file_1> <input_file_2>")
        print("Usage: python divide_bedgraph_values.py <input_file_1> <input_file_2> remove_zero_by_zero")
        print("Error: Please provide exactly two input files.")
        sys.exit(1)

    # check if the user wants to remove zero by zero
    if len(sys.argv) > 3 and sys.argv[3] == "remove_zero_by_zero":
        main(sys.argv[1], sys.argv[2], "divide_remove_0_by_0")
    else:
        main(sys.argv[1], sys.argv[2], "divide")
