import pandas as pd
import sys

# Written with Google Gemini

def process_tsv(file_path):
    """
    Reads a TSV file, separates it into two DataFrames based on 'keep' columns,
    and prints basic information about each.

    Args:
        file_path (str): The path to the TSV file.

    Returns:
        tuple: A tuple containing two pandas DataFrames:
               - df_all_keep_true: DataFrame where all 'keep' columns are True.
               - df_at_least_one_keep_false: DataFrame where at least one 'keep' column is False.
    """
    try:
        # (0) Read the file using pandas
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None, None

    # Identify 'keep' columns
    keep_cols = [col for col in df.columns if col.startswith("keep")]

    if not keep_cols:
        print("Warning: No columns starting with 'keep' found in the file.")
        return df, pd.DataFrame()  # Return original and an empty DataFrame

    # (1) Separate dataframe
    # DataFrame where all 'keep' columns are True
    df_all_keep_true = df[df[keep_cols].all(axis=1)]

    # DataFrame where at least one 'keep' column is False
    df_at_least_one_keep_false = df[~df[keep_cols].all(axis=1)]

    print("\n--- DataFrame where all 'keep' columns are True ---")
    print(f"Shape: {df_all_keep_true.shape}")
    if not df_all_keep_true.empty:
        print(f"First 5 rows:\n{df_all_keep_true.head()}")
        print("\nDescriptive statistics for numeric columns (gof, gofVerify, gofNoCut, gofNoCutVerify):\n"
              f"{df_all_keep_true[['gof', 'gofVerify', 'gofNoCut', 'gofNoCutVerify']].describe()}")
    else:
        print("DataFrame is empty.")

    print("\n--- DataFrame where at least one 'keep' column is False ---")
    print(f"Shape: {df_at_least_one_keep_false.shape}")
    if not df_at_least_one_keep_false.empty:
        print(f"First 5 rows:\n{df_at_least_one_keep_false.head()}")
        print("\nDescriptive statistics for numeric columns (gof, gofVerify, gofNoCut, gofNoCutVerify):\n"
              f"{df_at_least_one_keep_false[['gof', 'gofVerify', 'gofNoCut', 'gofNoCutVerify']].describe()}")
    else:
        print("DataFrame is empty.")

    return df_all_keep_true, df_at_least_one_keep_false

if __name__ == "__main__":

    # Our goal here is to compare a fit produced during our initial model fitting procedure
    # and verification of that fit later on.
    # If you don't know what this means, please ignore this script - it is part of a checking
    # step and is not relevant to the main analysis.

    if len(sys.argv) != 2:
        print("Usage: python gof_gofVerify_compare.py <file_path_tsv>")
        sys.exit(1)
    else:
        file_path = sys.argv[1]
        df_true, df_false = process_tsv(file_path)

        # You can now work with df_true and df_false separately
        # For example:
        # if df_true is not None:
        #     print("\nFurther processing on df_true...")
        #     # Your code here

        # if df_false is not None:
        #     print("\nFurther processing on df_false...")
        #     # Your code here