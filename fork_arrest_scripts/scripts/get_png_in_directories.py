import os
import sys
import pandas as pd


# Function to get all .png files from a folder
def get_png_files(folder):
    try:
        return [file for file in os.listdir(folder) if file.endswith('.png')]
    except FileNotFoundError:
        print(f"Folder {folder} not found.")
        return []


if __name__ == "__main__":

    # ensure data is piped in
    if sys.stdin.isatty():
        print("Goal: get all .png files from each folder and add a new row for each .png file")
        print("Usage: python get_png_in_directories.py < data.tsv > new_data.tsv")
        print("data.tsv should be a tab-separated file with a column named 'folder' and can have other columns")
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # Using pandas to read the csv data from stdin
    df = pd.read_csv(sys.stdin, delimiter='\t', comment='#')

    # Creating a new list to hold rows
    new_rows = []

    # get all .png files from each folder and add a new row for each .png file
    for _, row in df.iterrows():
        png_files = get_png_files(row['folder'])
        if not png_files:
            raise ValueError(f"No .png files found in {row['folder']}")
        else:
            for png_file in png_files:
                new_row = row.copy()
                new_row['plot_file'] = png_file
                new_rows.append(new_row)

    # Convert new rows to DataFrame and then write it to stdout
    new_df = pd.DataFrame(new_rows)
    new_df.to_csv(sys.stdout, sep='\t', index=False)
