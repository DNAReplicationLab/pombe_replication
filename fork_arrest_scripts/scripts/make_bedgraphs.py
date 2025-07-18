import sys
from user_input_tools import receive_pandas_and_params_from_user

if __name__ == "__main__":

    desc = ("""    Make bedgraphs given 1 + 3 or 1 + 4 column data
    separated by spaces or tabs. Lines starting with # are
    ignored. Headers are required (see below).

    Sample input:
    >> echo 'index contig start end val
    ind1 chrI 2 3 0.2
    ind1 chrI 1 2 0.1
    ind2 chrII 100 200 0.4
    ind2 chrII 200 300 0.5
    ind3 chrIV 1000 2000 0.4
    ind3 chrIV 1500 2500 0.5' | python make_bedgraphs.py --dir /tmp --ncol 4

    Sample output:
    >> cat /tmp/ind1.bedGraph
    track type=bedGraph name="ind1" description="ind1" visibility=full color=200,100,0 altColor=0,100,200 priority=20"""
            """ viewLimits=0.0:1.0
    chrI 1 2 0.1
    chrI 2 3 0.2
    >> cat /tmp/ind2.bedGraph
    track type=bedGraph name="ind2" description="ind2" visibility=full color=200,100,0 altColor=0,100,200 priority=20"""
            """ viewLimits=0.0:1.0
    chrII 100 200 0.4
    chrII 200 300 0.5
    >> cat /tmp/ind3.bedGraph
    cat: /tmp/ind3.bedGraph: No such file or directory

    Input: If 1 + 4 columns are supplied, the program interprets them
    as index,contig,start,end,numerical value. If 1 + 3 cols, interpretation
    is index,contig,start,numerical value with end = start + 1.
    Column headers must be index contig start val or index contig
    start end val. Columns can be in any order.

    Output: One unique file per index called whateverTheIndexIs.bedGraph
    in the given directory.

    NOTE: Numerical value must be b/w 0 and 1.

    NOTE: We check for invalid inputs like overlapping windows.
    If they are found, the corresponding index is ignored.
    The checks can be disabled using an input flag.
    
    NOTE: contig can be omitted if implicitContigFormat is used.
    But, column numbers still have to be set to 3 or 4.
    """)

    # get options and data
    df, args = receive_pandas_and_params_from_user([('dir', str, None, 'directory to put the bedgraphs'),
                                                    ('ncol', int, None, 'set to 3 or 4'),
                                                    ('comment', str, '', '(optional) a comment line in the header'),
                                                    ('suffix', str, '', '(optional) a suffix such that filenames are '
                                                                        'whateverTheIndexIsSuffix.bedGraph'),
                                                    ('implicitContigFormat', str, '',
                                                     '(optional) if index has delimiters and contains contig, '
                                                     'then contig can be inferred. e.g. if the index is '
                                                     'readID_contig_start_end, then set parameter to _2 to say '
                                                     'contig is at position 2 and delimiter is _. Overrides the '
                                                     'user-supplied contig column if present.'),
                                                    ('color', str, '200,100,0',
                                                     '(optional) 3 nums b/w 0 255 sep by ",". default 200,100,0'),
                                                    ],
                                                   [('noChecks', False, ('within an index, if some simple checks fail '
                                                                         '(e.g. start < end, no window overlap etc.) '
                                                                         'then that index is ignored unless this '
                                                                         'option is used. Do not use the option unless '
                                                                         'you are an experienced user. '))],
                                                   desc, sys.stdin,
                                                   required_columns=['index', 'start', 'val'],
                                                   additional_columns_allowed=True)

    df['valNotBw0and1'] = (df['val'] < 0) | (df['val'] > 1)

    # reject if numerical vals are not b/w 0 and 1
    if df.valNotBw0and1.any():
        raise ValueError('Only vals b/w 0 and 1!')

    # if three cols, make a fourth col
    if not (args.ncol == 3 or args.ncol == 4):
        raise ValueError("ncol must be 3 or 4")
    elif args.ncol == 3:
        df['end'] = df['start'] + 1

    for name, group_temp in df.groupby('index'):

        # copy the subset
        group = group_temp.copy()

        # sort the windows
        group.sort_values(by=['start'], inplace=True)

        if not args.noChecks:

            # do some checks
            group['prevEnd'] = group['end'].shift(1)
            group['stGeEn'] = group['start'] >= group['end']
            group['prevEnGSt'] = group['prevEnd'] > group['start']

            # go to the next index if failure
            if group.stGeEn.any() or group.prevEnGSt.any():
                continue

            # infer contig if option has been set
            if len(args.implicitContigFormat) > 0:
                group['contig'] = name.split(args.implicitContigFormat[0])[int(args.implicitContigFormat[1:]) - 1]

        with open(f"{args.dir}/{name}{args.suffix}.bedGraph", "w") as fp:

            # print bedgraph
            header = (f"track type=bedGraph name=\"{name}\""
                      f" description=\"{name}\""
                      f" visibility=full color={args.color} altColor=0,100,200"
                      " priority=20 viewLimits=0.0:1.0")

            if args.comment:
                print("#" + args.comment, file=fp)

            print(header, file=fp)
            print(group.to_csv(index=False, header=False,
                               columns=['contig', 'start', 'end', 'val'], sep=" "), file=fp)
