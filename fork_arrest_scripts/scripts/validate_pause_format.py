import sys
import pandas as pd
from DNAscentTools.modBAM_tools_additional import process_fork_index


def check_if_valid_pause_dataframe(df: pd.DataFrame) -> bool:
    """ Check that the input data frame is a valid pause dataframe.

    Args:
        df: pandas data frame to check

    Returns:
        True if the data frame is a valid pause dataframe, else False
    """

    # check that columns detectIndex and pauseSite are present
    if not all(col in df.columns for col in ['detectIndex', 'pauseSite']):
        return False

    # process the detectIndex column using process_fork_index
    try:
        df[['__read_id', '__contig', '__start', '__end', '__orn', '__dirn', '__startFork', '__endFork']] \
            = df['detectIndex'].apply(process_fork_index).tolist()
    except AssertionError:
        return False

    # set __score to 1000 if the entries under all columns that start with "keep" are True, else set to 0
    df['__score'] = 1000
    for col in df.columns:
        if col.startswith('keep'):
            df['__score'] = df['__score'] * df[col]

    # Check that pauseSite is within __startFork and __endFork for rows with __score equal to 1000
    if not all(row['__startFork'] <= row['pauseSite'] <= row['__endFork'] for _, row in
               df[df['__score'] == 1000].iterrows()):
        return False

    # if the 'pauseDuration' column is present, check that it is a positive integer for rows with __score equal to 1000
    if 'pauseDuration' in df.columns:
        if not all(row['pauseDuration'] > 0 for _, row in df[df['__score'] == 1000].iterrows()):
            return False

    return True


if __name__ == '__main__':

    # print usage if no input is provided
    if sys.stdin.isatty():
        print('Usage: cat <file> | python validate_pause_format.py')
        print('Usage: cat <file> | python validate_pause_format.py --required-columns <column1>,<column2>,...')
        sys.exit(1)

    # load pandas data frame from stdin, which is tab-separated columns with a header and comments start with '#'
    df_data = pd.read_csv(sys.stdin, sep='\t', comment='#', header="infer")

    # if --required-columns is provided, check that all required columns are present
    if len(sys.argv) > 1 and sys.argv[1] == '--required-columns':
        required_columns = sys.argv[2].split(',')
        if not all(col in df_data.columns for col in required_columns):
            print('invalid')
            sys.exit(1)

    if check_if_valid_pause_dataframe(df_data):
        print('valid')
    else:
        print('invalid')
