import sys
import os
from user_input_tools import receive_pandas_and_params_from_user
from DNAscentTools.modBAM_tools_additional import process_fork_index
from validate_pause_format import check_if_valid_pause_dataframe

if __name__ == "__main__":
    desc = """    Convert file containing results of pause site detection to bed file. 

    Input: Pipe in pause file contents (this is a plain text file).
           This file is tab-separated, comments start with '#', and the first non-comment line is column names.
           Only the following columns are processed here: detectIndex, and (optionally) keep*.
           - detectIndex is a string in the format readID_contig_start_end_orn_dirn_startFork_endFork.
               - Names are self-explanatory.
               - start, end, startFork, endFork are integers and represent positions on the reference genome.
               - startFork < endFork and start < end irrespective of fork and read direction.
               - orn is fwd/rev and dirn is L/R.
               - contig can contain underscores.
           - keep* columns are optional and are boolean. If they are present and upon user request,
             they are used to filter the data. If any of the keep* columns is False, the pause is discarded.

    Output: Space-separated data in eight column format with no header or column names to two files in the output
            directory. The columns are contig, startFork, endFork, readID, contig, start, end, orn.
            Left forks go to leftForks_DNAscent_forkSense.bed and right forks go to rightForks_DNAscent_forkSense.bed.
            NOTE: these are not true bed files, we are following the notation of DNAscent.
            NOTE: Depending on the input data, these may not even be DNAscent forkSense calls.
                  They are just the format DNAscent uses.

    Sample usage:
        < pauseFile.txt python convert_pause_file_to_forkSense_calls.py --outputDir outputDir
        < pauseFile.txt python convert_pause_file_to_forkSense_calls.py --outputDir outputDir --discardUnkeptPauses
    """

    # get options and data
    df, args = receive_pandas_and_params_from_user(
        [('outputDir', str, None, 'output directory to send forkSense calls to')],
        [('discardUnkeptPauses', False, 'discard pauses where any keep* column is False')],
        desc, sys.stdin,
        required_columns=['detectIndex'],
        additional_columns_allowed=True,
        force_tab_separator=True)

    # check that the output directory exists
    if not os.path.exists(args.outputDir):
        raise Exception('Output directory does not exist: ' + args.outputDir)

    # check that the input data is a valid pause dataframe
    if not check_if_valid_pause_dataframe(df):
        raise Exception('Input data is not a valid pause file')

    # process the detectIndex column using process_fork_index
    df[['__read_id', '__contig', '__start', '__end', '__orn', '__dirn', '__startFork', '__endFork']] \
        = df['detectIndex'].apply(process_fork_index).tolist()

    # set __score to 1000 if the entries under all columns that start with "keep" are True, else set to 0
    df['__score'] = 1000
    for col in df.columns:
        if col.startswith('keep'):
            df['__score'] = df['__score'] * df[col]

    if args.discardUnkeptPauses:
        df = df[df['__score'] == 1000]

    # print the following columns in the following order to outputDirectory/leftForks_DNAscent_forkSense.bed
    # for all entries with __dirn == 'L' and to outputDirectory/rightForks_DNAscent_forkSense.bed for all entries
    # with __dirn == 'R'
    output_cols = ['__contig', '__startFork', '__endFork', '__read_id', '__contig', '__start', '__end', '__orn']
    df[df['__dirn'] == 'L'][output_cols].to_csv(args.outputDir + '/leftForks_DNAscent_forkSense.bed', sep=' ',
                                                header=False, index=False)
    df[df['__dirn'] == 'R'][output_cols].to_csv(args.outputDir + '/rightForks_DNAscent_forkSense.bed', sep=' ',
                                                header=False, index=False)
