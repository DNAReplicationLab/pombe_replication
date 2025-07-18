import sys
import argparse
import pandas as pd
from collections.abc import Callable
from typing import IO
from io import StringIO


def receive_pandas_and_params_from_user(param_list: list[tuple[str, Callable[[str], any], any, str]],
                                        flag_list: list[tuple[str, bool, str]],
                                        desc: str, fp: IO[str], required_columns: list[str],
                                        additional_columns_allowed: bool = True,
                                        allow_blank_piped_input: bool = False,
                                        force_tab_separator: bool = False) \
        -> tuple[pd.DataFrame, argparse.Namespace]:
    """ Gets data piped in by user and command line options

    Args:
        param_list: each entry tuple of parameter name, parameter type,
            default value (if no default and must be given by user, set to None), description
        flag_list: optional flags. each entry is a tuple of flag name, default value, description. if flag is set
            by user, then value of flag = (not default). flag can only store boolean.
        desc: description of program
        fp: file-like object, set to sys.stdin unless you are an advanced user
        required_columns: columns required in dataframe
        additional_columns_allowed: (default T) allow columns other than the ones required
        allow_blank_piped_input: (default F) allow program to be run with no piped input
        force_tab_separator: (default F) force tab separator instead of (generic) whitespace

    Returns:
        tuple of pandas dataframe, argparse namespace

    """

    if required_columns is None:
        required_columns = []

    # get options if param_list and/or flag_list are given
    if param_list == [] and flag_list == []:
        args = argparse.Namespace()
    else:
        parser = argparse.ArgumentParser(description=desc,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)

        for k in param_list:
            parser.add_argument(f'--{k[0]}', type=k[1], required=(k[2] is None), default=k[2], help=k[3])

        for k in flag_list:
            def get_action(x): return 'store_false' if x else 'store_true'
            parser.add_argument(f'--{k[0]}', required=False, default=k[1], action=get_action(k[1]), help=k[2])

        args = parser.parse_args()

    # read piped in data if available, else check if a blank input is fine, else fail.
    if not fp.isatty():
        inp_text = fp.read()
    elif allow_blank_piped_input:
        return pd.DataFrame(), args
    else:
        raise NotImplementedError("Please pipe in inputs")

    # get pandas dataframe
    if force_tab_separator:
        df_val = pd.read_csv(StringIO(inp_text), sep='\t',
                             comment="#", index_col=None)
    else:
        df_val = pd.read_csv(StringIO(inp_text), sep='\s+',
                             comment="#", index_col=None)

    # check if required columns are present
    if not all(k in df_val.columns for k in required_columns):
        raise ValueError("Required columns are missing")

    # disallow additional columns if necessary
    if not (additional_columns_allowed or len(df_val.columns) == len(required_columns)):
        raise ValueError("Additional columns not allowed")

    # return data
    return df_val, args


if __name__ == "__main__":
    description = """    Test if getting user input using pipes and command line works.
    * Pipe in one column of numbers with header a and use --num on command line and specify
      another number.
    * Observe if you receive an output where the column has been multiplied
      with the number.
    """

    df, args_received = receive_pandas_and_params_from_user([("num", int, None, "integer input")], [],
                                                            description, sys.stdin, ['a'])

    df['a'] = df['a'] * args_received.num

    print(df.to_csv(index=False, header=True), end="")
