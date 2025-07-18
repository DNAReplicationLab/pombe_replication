#!/usr/bin/env python

import json
import sys


def convert_json_to_latex(input_json: str) -> str:
    """ Converts a json object containing array of values to a latex list

    Args:
        input_json: a string like [{"value": 1, "desc": "", "latex_desc": "" },...]

    Returns:
        a string like "\begin{itemize}
                       \item [desc] latex_desc \\
                       value
                       ...
                       \end{itemize}"
                        
    """
    data = json.loads(input_json)

    latex_snippet = "\\begin{itemize}\n"

    for item in data:
        if all(key in item for key in ['desc', 'value', 'latex_desc']):

            # replace '_' with '\_' to avoid subscript in item['desc']
            item['desc'] = item['desc'].replace('_', '\_')

            latex_snippet += '\\item{}' + '[{}] {} \\\\\n'.format(item['desc'], item['latex_desc'])
            latex_snippet += str(item['value']) + '\n'

    latex_snippet += "\\end{itemize}"

    return latex_snippet


if __name__ == "__main__":

    # read piped in data
    if sys.stdin.isatty():
        raise NotImplementedError("Pipe in inputs please.")
    input_json = sys.stdin.read()

    latex_str = convert_json_to_latex(input_json)
    print(latex_str)
