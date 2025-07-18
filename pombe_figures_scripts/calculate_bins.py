import pandas as pd
import sys

def calculate_average_bins(df, bin_size=1000):
    result_bins = []

    for chromosome, data in df.groupby('chromosome'):
        bins = []
        start_bin = 0  # Initialize start_bin at 0
        sum_count = 0
        count_in_bin = 0

        # Determine the maximum position for this chromosome
        max_position = data['position'].max()

        while start_bin <= max_position:
            # Calculate end_bin based on start_bin and bin_size
            end_bin = start_bin + bin_size 

            # Filter data within the current bin range
            bin_data = data[(data['position'] >= start_bin) & (data['position'] <= end_bin)]

            if not bin_data.empty:
                # Calculate average count for the current bin
                average_count = bin_data['count'].mean()
                rounded_count = round(average_count)  # Round to nearest integer

                # Append bin details to bins list
                bins.append({
                    'chromosome': chromosome,
                    'start_bin': start_bin,
                    'end_bin': end_bin,
                    'average_count': average_count,
                    'rounded_count': rounded_count,
                })

            # Move to the next bin
            start_bin += bin_size

        result_bins.extend(bins)

    return result_bins

def main(input_file, output_file):
    try:
        df = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'position', 'count'])
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{input_file}' not found.")

    # Calculate average bins with default bin_size=1000
    result_bins = calculate_average_bins(df, bin_size=1000)

    # Convert list of dictionaries to DataFrame
    result_df = pd.DataFrame(result_bins)

    # Write result DataFrame to output file
    result_df.to_csv(output_file, sep='\t', index=False)

    print(f"Processing complete. Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: ./calculate_bins.py <input_file> <output_file>")

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    main(input_file, output_file)
