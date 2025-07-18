import sys
import pandas as pd
import re


def check_if_valid_summary_dataframe(df: pd.DataFrame) -> bool:
    """ Check that the input data frame is a valid summary dataframe.

    Args:
        df: pandas data frame to check

    Returns:
        True if the data frame is a valid pause dataframe, else False
    """

    # check that a list of required columns are present
    col_list = ["expected", "sd", "observed", "coverage", "coverage_sd", "division", "dataset",
                "feature", "analysis_label", "relative_direction"]
    if not all(col in df.columns for col in col_list):
        print("missing required columns")
        return False

    # check that the expected, sd, observed, coverage, and coverage_sd columns only contain numeric values
    if not all(df[col].apply(lambda x: isinstance(x, (int, float))).all() for col in
               ["expected", "sd", "observed", "coverage", "coverage_sd"]):
        print("expected, sd, observed, coverage, and coverage_sd columns must only contain numeric values")
        return False

    # check that the relative_direction column only contains the values 'all', 'co-directional', and 'head-on'
    if not all(df['relative_direction'].isin(['all', 'co-directional', 'head-on'])):
        print("relative_direction column must only contain the values 'all', 'co-directional', and 'head-on'")
        return False

    # check that the division column is either "NA" or "VS_integer" or "LS_integer" or "VS_integer_LS_integer"
    # using regular expressions
    if not df['division'].apply(lambda x: (isinstance(x, str) and
                                           (x == "NA" or re.match(r'VS_\d+', x) or
                                            re.match(r'LS_\d+', x) or
                                            re.match(r'VS_\d+_LS_\d+', x)))
                                or pd.isna(x)).all():
        print("division column must be either 'NA' or 'VS_integer' or 'LS_integer' or 'VS_integer_LS_integer'")
        return False

    return True


if __name__ == '__main__':

    # print usage if no input is provided
    if sys.stdin.isatty():
        print('Usage: cat <file> | python validate_summary_file.py')
        sys.exit(1)

    # load pandas data frame from stdin, which is tab-separated columns with a header and comments start with '#'
    df_data = pd.read_csv(sys.stdin, sep='\t', comment='#', header="infer")

    if check_if_valid_summary_dataframe(df_data):
        print('valid')
    else:
        print('invalid')
