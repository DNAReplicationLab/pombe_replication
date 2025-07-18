import argparse
import sys
import pandas as pd
from DNAscentTools.modBAM_tools_additional import process_fork_index


def main(data: pd.DataFrame, forkLocThres: int, pauseLocThres: int) -> pd.DataFrame:
    """ Mark the forks and pauses that are close to the ends of the alignment.

    Args:
        data: A pandas dataframe obtained from the input pause file. Must contain a detectIndex and a pauseSite column.
        forkLocThres: The threshold in base pairs for the distance between the fork and the ends.
        pauseLocThres: The threshold in base pairs for the distance between the pause site and the ends.

    Returns:
        A pandas dataframe with two additional columns: keep_if_fork_not_close_to_end and keep_if_pause_not_close_to_end
    """

    # Check that required columns exist
    required_columns = ['detectIndex', 'pauseSite']
    for col in required_columns:
        if col not in data.columns:
            raise ValueError(f"Required column '{col}' not found in the input file.")

    df[['__read_id', '__contig', '__start', '__end', '__orn', '__dirn', '__startFork', '__endFork']] \
        = df['detectIndex'].apply(process_fork_index).tolist()

    # Check if optional columns exist, else create them
    optional_columns = ['keep_if_fork_not_close_to_end', 'keep_if_pause_not_close_to_end']
    for col in optional_columns:
        if col not in data.columns:
            data[col] = None

    # Verify the integrity of the input data
    if not data.apply(lambda row: row['__start'] <= row['__startFork'] < row['__endFork'] <= row['__end'], axis=1).all():
        raise ValueError("All entries must satisfy start <= startFork < endFork <= end")

    if not data.apply(lambda row: row['__startFork'] < row['pauseSite'] < row['__endFork'], axis=1).all():
        raise ValueError("All entries must satisfy startFork < pauseSite < endFork")

    # Process the data
    data['keep_if_fork_not_close_to_end'] = data.apply(
        lambda row: not (row['__startFork'] - row['__start'] <= forkLocThres or
                         row['__end'] - row['__endFork'] <= forkLocThres),
        axis=1
    )

    data['keep_if_pause_not_close_to_end'] = data.apply(
        lambda row: not (row['pauseSite'] - row['__start'] <= pauseLocThres or
                         row['__end'] - row['pauseSite'] <= pauseLocThres),
        axis=1
    )

    # drop all columns starting with '__'
    data = data.loc[:, ~data.columns.str.startswith('__')]

    # return the processed data
    return data


if __name__ == "__main__":
    desc = r"""
        Sample usage shown below:

        ```bash
        python mark_forks_pauses_at_align_ends.py pauseFile > pauseFile_2
        python mark_forks_pauses_at_align_ends.py --forkLocThresKb 10 pauseFile > pauseFile_2
        python mark_forks_pauses_at_align_ends.py --pauseLocThresKb 10 pauseFile > pauseFile_2
        python mark_forks_pauses_at_align_ends.py --forkLocThresKb 10 --pauseLocThresKb 10 pauseFile > pauseFile_2
        ```

        - pauseFile is tab-separated with column names, and comments that start with '#'.
        - Check that columns detectIndex, pauseSite exist.
        - Process detectIndex to extract start, end, startFork, endFork.
        - Check if columns keep_if_fork_not_close_to_end, keep_if_pause_not_close_to_end exist, else create them.
        - all entries satisfy start <= startFork < endFork <= end
        - all entries satisfy startFork < pauseSite < endFork
        - if startFork - start or end - endFork are within forkLocThresKb * 1000, then set the corresponding entr in 
          column keep_if_fork_not_close_to_end to False, else set it to True.
        - if pauseSite - start or end - pauseSite are within pauseLocThresKb * 1000, then the corresponding entry in 
          column keep_if_pause_not_close_to_end to False, else set it to True.

        Write output to stdout, which is the same as the input file, but with two additional columns as described above.

    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("pauseFile", type=str, help="Path to the input file with tab-separated data.")
    parser.add_argument("--forkLocThresKb", type=int, default=0,
                        help="Threshold in kilobases for the distance between the fork and the ends.")
    parser.add_argument("--pauseLocThresKb", type=int, default=0,
                        help="Threshold in kilobases for the distance between the pause site and the ends.")
    args = parser.parse_args()

    # Read the input file using pandas
    df = pd.read_csv(args.pauseFile, sep='\t', comment='#')

    main(df, args.forkLocThresKb * 1000, args.pauseLocThresKb * 1000).to_csv(sys.stdout, sep='\t', index=False)
