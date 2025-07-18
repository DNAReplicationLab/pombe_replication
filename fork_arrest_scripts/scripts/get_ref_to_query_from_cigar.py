import sys
import pandas as pd
from io import StringIO
from DNAscentTools.modBAM_tools_additional import cigar_to_ref_to_query_tbl

if __name__ == "__main__":

    desc = """    Obtain reference to query coordinate mapping.

    Input: Pipe in three columns w headers index, posOnRefStart,
    cigarStr in any order. Data separator should be spaces or tabs.
    index is any index that uniquely marks each line,
    posOnRefStart is an integer coordinate for start of
    position on reference, cigarStr self-explanatory.
    
    Output is to stdout and 3 columns are present: index,
    refCoord, queryCoord. Data separator is spaces.
    Column headers are not output.
    
    NOTE: Really long CIGAR operations are continued in the CG: tag by minimap2
    when run with suitable options. This script does not handle those.
    So if you have a CG: tag in your BAM file, either you should not use this script,
    or find a way to incorporate the CG: tag into the input data fed to this script.
    For more details, look at the section 'Working with >65535 CIGAR operations'
    at https://github.com/lh3/minimap2 .
    """
    # read piped in data
    if not sys.stdin.isatty():
        inpText = sys.stdin.read()
    else:
        raise NotImplementedError("Pipe in inputs please.")

    # load data
    dfVal = pd.read_csv(StringIO(inpText), sep="\s+",
                        comment="#", index_col='index')

    for index, k in dfVal.iterrows():
        for m in cigar_to_ref_to_query_tbl(k.cigarStr, k.posOnRefStart):
            print(f"{index} {m[0]} {m[1]}")
