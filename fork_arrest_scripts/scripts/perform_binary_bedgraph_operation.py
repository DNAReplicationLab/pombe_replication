import sys
import os
import pandas as pd


def perform_operation(df_join: pd.DataFrame, operation: str) -> pd.DataFrame:
    """ Perform operation on corresponding values from two bedgraph files.

    Args:
        df_join: pandas.DataFrame with columns "val_a" and "val_b"
        operation: string specifying operation: "add", "subtract", "multiply", "divide" or "divide_remove_0_by_0".
                   divide_remove_0_by_0 is a special case where any 0/0 is removed.

    Returns:
        pandas.DataFrame with a result column added
    """
    if operation == "add":
        df_join["result"] = df_join["val_a"] + df_join["val_b"]
    elif operation == "subtract":
        df_join["result"] = df_join["val_a"] - df_join["val_b"]
    elif operation == "multiply":
        df_join["result"] = df_join["val_a"] * df_join["val_b"]
    elif operation == "divide":
        if (df_join["val_b"] == 0).any():
            print("Error: Division by zero encountered in val_b.")
            sys.exit(1)
        df_join["result"] = df_join["val_a"] / df_join["val_b"]
    elif operation == "divide_remove_0_by_0":
        # check for non-zero/zero
        if df_join.apply(lambda x: x['val_b'] == 0 & x['val_a'] != 0, axis=1).any():
            print("Error: Non-zero divided by zero encountered.")
            sys.exit(1)
        df_join["result"] = df_join["val_a"] / df_join["val_b"]
        df_join.dropna(inplace=True)
    else:
        print(f"Error: Invalid operation {operation}.")
        sys.exit(1)

    return df_join


def process_bedgraphs(file_a: str, file_b: str, operation: str) -> pd.DataFrame:
    """ Process two bedgraph files, perform a specified operation on the values, and return the results.

    Args:
        file_a: path to bedgraph file
        file_b: path to bedgraph file
        operation: math operation to perform on the values, can be "add", "subtract", "multiply", or "divide"

    Returns:
        pandas.DataFrame with columns "contig", "start", "end", and "result"
    """
    # Read input files
    df_a = pd.read_csv(file_a, sep=" ", names=["contig", "start", "end", "val"], comment="#")
    df_b = pd.read_csv(file_b, sep=" ", names=["contig", "start", "end", "val"], comment="#")

    # Join data frames
    df_join = df_a.merge(df_b, on=["contig", "start", "end"], suffixes=("_a", "_b"))

    # raise error if any start > end
    if (df_join["start"] >= df_join["end"]).any():
        print("Error: start is greater than or equal to end in some rows.")
        sys.exit(1)

    # Perform specified operation
    df_result = perform_operation(df_join, operation)

    # Return the result data frame
    return df_result[["contig", "start", "end", "result"]]


def main(file_1: str, file_2: str, operation: str) -> None:
    """ Process two bedgraph files, perform a specified operation on the values, and output the results to stdout.

    Args:
        file_1: path to bedgraph file
        file_2: path to bedgraph file
        operation: math operation to perform, can be "add", "subtract", "multiply", "divide" or "divide_remove_0_by_0"

    Returns:
        None
    """
    # Check that input files exist
    for file in [file_1, file_2]:
        if not os.path.isfile(file):
            print(f"Error: File {file} does not exist.")
            sys.exit(1)

    # Process input files
    print(f"# {operation} {file_1} and {file_2}")
    df_result = process_bedgraphs(file_1, file_2, operation)

    # Output results
    df_result[["contig", "start", "end", "result"]].to_csv(sys.stdout, sep=" ", header=False, index=False)


if __name__ == "__main__":

    # check inputs
    if len(sys.argv) != 4:
        print("Error: Please provide exactly two input files and one operation.")
        print("Usage: python perform_binary_bedgraph_operation.py file_1 file_2 operation")
        print("Output is written to stdout.")
        print("Operation can be add, subtract, multiply, or divide.")
        print("Example: python perform_binary_bedgraph_operation.py file_1.bedgraph file_2.bedgraph add")
        sys.exit(1)

    # run main logic function
    main(sys.argv[1], sys.argv[2], sys.argv[3])
