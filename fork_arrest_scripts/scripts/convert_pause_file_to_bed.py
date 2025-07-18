from user_input_tools import receive_pandas_and_params_from_user
from DNAscentTools.modBAM_tools_additional import process_fork_index
from validate_bed_format import validate_first_three_columns, validate_next_three_columns
import sys

if __name__ == "__main__":

    desc = """    Convert file containing results of pause site detection to bed file. 

    Input: Pipe in pause file contents (this is a plain text file).
           This file is tab-separated, comments start with '#', and the first non-comment line is column names.
           Only the following columns are processed here: detectIndex, pauseSite, keep*, pauseDuration, paramStrNoPause.
           - First two columns are required, "keep*" columns are optional, and pauseDuration and paramStrNoPause are
             processed only if the appropriate options are used.
           - detectIndex is a string in the format readID_contig_start_end_orn_dirn_startFork_endFork.
               - Names are self-explanatory.
               - start, end, startFork, endFork are integers and represent positions on the reference genome.
               - startFork < endFork and start < end irrespective of fork and read direction.
               - orn is fwd/rev and dirn is L/R.
               - contig can contain underscores.
           - pauseSite is an integer and represents a position on the reference genome.
           - keep* are boolean values "True" or "False" and represent whether the pause site should be kept or not
             depending on the criteria used. (* is a wildcard and can be any string, so columns like keep_ab,
             keep_cd etc. fall under this category. NOTE that ab cd in this example are fictional).
           - If keep* columns are not present, then all pause sites are kept.
           - If outputSgmStartToCenter is used, then the paramNoStrPause column is needed.
             This column has the format high_low_offset_width where all four are numbers.
                high and low are the high and low levels of the sigmoid, offset is the
                center point of the sigmoid along the reference genome, and width is the
                sigmoid width with the sign convention of + for right-moving fork and - for left.
           - If outputDurationAsScore is used, then the pauseDuration column is needed.
             This is the duration of the pause, usually in kb, but this script is unit-agnostic.
             If a pause is marked to be discarded, then its duration may be undefined.
             So, the pauseDuration column is output as score only if the discardUnkeptPauses option is used.

    Output: Tab-separated data in BED6+N format with comments starting with '#' are output to stdout.
            The standard six BED columns are: chrom, chromStart, chromEnd, name, score, strand.
            The other +N columns are whatever else is in the input file, in the same order.
            Score is set to 1000 or 0 depending on whether the pause site is kept or not respectively, or
            set to pauseDuration if outputDurationAsScore is used.
            chromStart and chromEnd can be pauseSite and pauseSite+1 respectively, or
            other columns like startFork and endFork depending on user options.

    Sample usage:
        < pauseFile.txt python convert_pause_file_to_bed.py > pauseFile.bed
        < pauseFile.txt python convert_pause_file_to_bed.py --outputForks > pauseFile.bed 
        < pauseFile.txt python convert_pause_file_to_bed.py --outputForks --LRtoPlusMinus > pauseFile.bed 
    """

    # get options and data
    df, args = receive_pandas_and_params_from_user([],
                                                   [('outputForks', False, ('output startFork and endFork instead of '
                                                                            'pauseSite and pauseSite + 1.')),
                                                    ('LRtoPlusMinus', False, ('use fork direction to mark '
                                                                              'strand as + (R) or - (L)')),
                                                    ('outputAlignments', False, ('output start and end instead of '
                                                                                 'pauseSite and pauseSite + 1.')),
                                                    ('LeadLagToPlusMinus', False, ('mark strand as - (lagging) or + ('
                                                                                   'leading)')),
                                                    ('discardUnkeptPauses', False, ('discard pauses where any keep* '
                                                                                    'column is False')),
                                                    ('outputBed6Plus1WherePlus1isAlignStrand', False,
                                                     'output 6+1 columns where +1 is alignment strand (+/-)'),
                                                    ('outputSgmStartToCenter', False,
                                                     'output start of sigmoid till its center instead of pauseSite and '
                                                     'pauseSite + 1. Need paramNoStrPause column, see comments in the '
                                                     'script for more details.'),
                                                    ('outputDurationAsScore', False,
                                                     'output duration of pause as score instead of 1000/0. Must use '
                                                     'with --discardUnkeptPauses'),
                                                    ('outputNormStepAsScore', False,
                                                     'output data from abs_step_normalized_by_sd of pause as score '
                                                     'instead of 1000/0. Must use with --discardUnkeptPauses. '
                                                     'This column is reported in pause files generated by the rDNA '
                                                     'pause pipeline.'),
                                                    ('outputDetectIndexAsFourthColumn', False,
                                                     'outputs the detect index instead of the read id '
                                                     'as the fourth column in the bed file.')
                                                    ],
                                                   desc, sys.stdin,
                                                   required_columns=['detectIndex', 'pauseSite'],
                                                   additional_columns_allowed=True,
                                                   force_tab_separator=True)

    # if the input file is empty, exit
    if len(df.index) == 0:
        sys.exit(0)

    columns = df.columns.tolist()

    # check that outputForks and outputAlignments are not both used
    if (args.outputForks and args.outputAlignments) or (args.outputForks and args.outputSgmStartToCenter) or \
            (args.outputAlignments and args.outputSgmStartToCenter):
        raise ValueError('Error: Use only one of --outputForks, --outputAlignments, --outputSgmStartToCenter')

    # check that LeadLagToPlusMinus is not used with LRtoPlusMinus
    if args.LeadLagToPlusMinus and args.LRtoPlusMinus:
        raise ValueError('Error: --LeadLagToPlusMinus and --LRtoPlusMinus cannot be used together.')

    # check that outputDurationAsScore and outputNormStepAsScore are not used together
    # or used without discardUnkeptPauses
    if (args.outputDurationAsScore or args.outputNormStepAsScore) and not args.discardUnkeptPauses:
        raise ValueError('Error: --outputDurationAsScore or --outputNormStepAsScore must be used with '
                         '--discardUnkeptPauses')
    if args.outputDurationAsScore and args.outputNormStepAsScore:
        raise ValueError('Error: --outputDurationAsScore and --outputNormStepAsScore cannot be used together')

    # if outputDurationAsScore is used, then pauseDuration column must be present
    if args.outputDurationAsScore and 'pauseDuration' not in columns:
        raise ValueError('Error: pauseDuration column not found in input file, see --outputDurationAsScore')

    # if outputNormStepAsScore is used, then abs_step_normalized_by_sd column must be present
    if args.outputNormStepAsScore and 'abs_step_normalized_by_sd' not in columns:
        raise ValueError('Error: abs_step_normalized_by_sd column not found in input file, see --outputNormStepAsScore')

    # set __score to 1000 if the entries under all columns that start with "keep" are True, else set to 0
    df['__score'] = 1000
    for col in df.columns:
        if col.startswith('keep'):
            df['__score'] = df['__score'] * df[col]

    # process the detectIndex column using process_fork_index
    df[['__read_id', '__contig', '__start', '__end', '__orn', '__dirn', '__startFork', '__endFork']] \
        = df['detectIndex'].apply(process_fork_index).tolist()

    # Check that pauseSite is within __startFork and __endFork for rows with __score equal to 1000, else raise error
    if not all(row['__startFork'] <= row['pauseSite'] <= row['__endFork'] for _, row in
               df[df['__score'] == 1000].iterrows()):
        raise ValueError('Error: pauseSite is not within startFork and endFork for rows with valid pauses')

    # convert __orn values to + or - depending on whether they are fwd or rev
    df['__orn'] = df['__orn'].apply(lambda x: '+' if x == 'fwd' else '-')
    align_strand = df['__orn'].tolist()

    # if --LRtoPlusMinus or --LeadLagToPlusMinus is used, then set __orn accordingly
    if args.LRtoPlusMinus:

        # set __orn to + if dirn is R, else set to -
        df['__orn'] = df['__dirn'].apply(lambda x: '+' if x == 'R' else '-')
        print('# strand of read (- or +) is set by fork direction (L or R)')

    elif args.LeadLagToPlusMinus:

        # FLAG: LEAD LAG DISTINCTION
        # set __orn to + if __dirn is R and __orn is fwd or __dirn is L and __orn is rev, else set to -
        df['__orn'] = df.apply(lambda x: '+' if ((x['__dirn'] == 'R' and x['__orn'] == '+') or
                                                 (x['__dirn'] == 'L' and x['__orn'] == '-')) else '-', axis=1)
        print('# strand of read is set by leading (+) or lagging (-) strand synthesis')
        print('# leading means fork direction, read direction = L,- or R,+')
        print('# lagging means fork direction, read direction = L,+ or R,-')

    # set start and end columns depending on whether --outputForks, --outputAlignments, or neither is used
    if args.outputForks:
        df['__startOutput'] = df['__startFork']
        df['__endOutput'] = df['__endFork']
        print('# start and end columns are startFork and endFork respectively')
    elif args.outputAlignments:
        df['__startOutput'] = df['__start']
        df['__endOutput'] = df['__end']
        print('# start and end columns are start and end of the alignment respectively')
    elif args.outputSgmStartToCenter:
        if "paramStrNoPause" not in columns:
            raise ValueError('Error: paramStrNoPause column not found in input file, see --outputSgmStartToCenter')

        # remove rows where paramStrNoPause has an invalid value
        df = df[df['paramStrNoPause'].apply(lambda x: isinstance(x, str))]

        # extract offset and width from paramStrNoPause column
        df['__offset'] = df['paramStrNoPause'].apply(lambda x: int(float(x.split('_')[2]))).tolist()
        width = df['paramStrNoPause'].apply(lambda x: int(float(x.split('_')[3]))).tolist()

        # check that extracted widths match with width column. this is strictly not necessary, but is a good check.
        width_check = df['width'].tolist()
        assert all([w == w_check for w, w_check in zip(width, width_check)])

        # set start and end in the bed file to sigmoid locations and clean up
        df['__startOutput'] = df[["__start", "__offset", "width"]].apply(lambda x: x[0] if x[2] > 0 else x[1], axis=1)
        df['__endOutput'] = df[["__offset", "__end", "width"]].apply(lambda x: x[0] if x[2] > 0 else x[1], axis=1)
        del df['__offset']

        # set __startOutput to __start if __startOutput is less than __start
        df['__startOutput'] = df[["__start", "__startOutput"]].apply(lambda x: x[0] if x[1] < x[0] else x[1], axis=1)

        # set __endOutput to __end if __endOutput is greater than __end
        df['__endOutput'] = df[["__end", "__endOutput"]].apply(lambda x: x[0] if x[1] > x[0] else x[1], axis=1)

        # remove invalid rows
        df = df[(0 <= df['__startOutput']) & (df['__startOutput'] <= df['__endOutput'])]
    else:
        df['__startOutput'] = df['pauseSite']
        df['__endOutput'] = df['pauseSite'] + 1
        print('# start and end columns are pauseSite and pauseSite + 1 respectively')

    if args.outputDurationAsScore:
        print('# pause duration is output as score and pauses marked for discard have been discarded')
    elif args.outputNormStepAsScore:
        print('# abs_step_normalized_by_sd is output as score and pauses marked for discard have been discarded')
    else:
        print('# score of 1000/0 means pause detected/pause discarded')
    print('# in our method, a fork can have zero or one pause sites')
    print('# format is BED6+N, where N is the number of additional columns in the input file or ')
    print('#     N = 1 where column is alignment strand, depending on user inputs to the script')

    # delete the __start, __end, __dirn columns
    del df['__start']
    del df['__end']
    del df['__dirn']
    del df['__startFork']
    del df['__endFork']

    # if outputDetectIndexAsFourthColumn is used, then set __read_id to detectIndex
    if args.outputDetectIndexAsFourthColumn:
        df['__read_id'] = df['detectIndex']

    # reorder columns
    canonical_bed_cols = ['__contig', '__startOutput', '__endOutput', '__read_id', '__score', '__orn']
    bed_cols = canonical_bed_cols.copy()
    bed_cols.extend([col for col in columns if col not in bed_cols])
    df = df.reindex(sorted(df.columns, key=lambda x: bed_cols.index(x) if x in bed_cols else len(bed_cols)), axis=1)

    # validate bed file before output
    if not all(validate_first_three_columns(k['__contig'], str(k['__startOutput']), str(k['__endOutput']))
               for k in df.to_dict('records')):
        raise ValueError('Error: Invalid values in first three columns')

    # pass __read_id, __score, __orn to validate_next_three_columns and ensure all values are True
    if not all(validate_next_three_columns(k['__read_id'], str(k['__score']), k['__orn'])
               for k in df.to_dict('records')):
        raise ValueError('Error: Invalid values in next three columns')

    if args.outputBed6Plus1WherePlus1isAlignStrand:
        df['alignStrand'] = align_strand

    if args.discardUnkeptPauses:
        df = df[df['__score'] == 1000]
        if args.outputDurationAsScore:
            df['__score'] = df['pauseDuration']
        if args.outputNormStepAsScore:
            df['__score'] = df['abs_step_normalized_by_sd']

    # output columns separated by tabs to stdout
    if not args.outputBed6Plus1WherePlus1isAlignStrand:
        df.to_csv(sys.stdout, sep='\t', index=False, header=False, na_rep='NA')
    else:
        df[canonical_bed_cols + ['alignStrand']].to_csv(sys.stdout, sep='\t', index=False,
                                                        header=False, na_rep='NA')
