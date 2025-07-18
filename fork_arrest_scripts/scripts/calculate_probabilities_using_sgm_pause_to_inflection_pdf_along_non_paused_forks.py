import sys
import os
from fork_arrest_tools import DNASequenceIterator, process_param_str
from DNAscentTools.modBAM_tools_additional import process_fork_index
from user_input_tools import receive_pandas_and_params_from_user
import scipy.stats as st
import numpy as np
import pandas as pd

if __name__ == "__main__":

    desc = """    Calculate probabilities of pause site detection along non-paused forks using stats from paused forks. 

    Input: Pipe in pause file contents (this is a plain text file).
           This file is tab-separated, comments start with '#', and the first non-comment line is column names.
           Only the following columns are processed here: detectIndex, pauseSite, paramStrNoPause, keep*.
           - First three columns are required, "keep*" columns are optional.
           - detectIndex is a string in the format readID_contig_start_end_orn_dirn_startFork_endFork.
               - Names are self-explanatory.
               - start, end, startFork, endFork are integers and represent positions on the reference genome.
               - startFork < endFork and start < end irrespective of fork and read direction.
               - orn is fwd/rev and dirn is L/R.
               - contig can contain underscores.
           - pauseSite is an integer and represents a position on the reference genome.
           - paramNoStrPause has the format high_low_offset_width where all four are numbers.
                high and low are the high and low levels of the sigmoid, offset is the
                center point of the sigmoid along the reference genome, and width is the
                sigmoid width with the sign convention of + for right-moving fork and - for left.
           - keep* are boolean values "True" or "False" and represent whether the pause site should be kept or not
             depending on the criteria used. (* is a wildcard and can be any string, so columns like keep_ab,
             keep_cd etc. fall under this category. NOTE that ab cd in this example are fictional).
           - If keep* columns are not present, then all pause sites are kept.
           
           Also set options on the command line.

    Output: Tab-separated data with header with column names and comments starting with '#' are output to stdout.
            The columns are detectIndex, pauseSite, probability, earlyPieceDurn latePieceDurn.
            - First column has already been described in the input section.
            - Second column is candidate pause sites along the fork.
            - Third column is the probability of finding pauses at the candidate pause site.
            - Fourth column is the distance of pause site from that end of the fork that has high BrdU.
            - Fifth column is the distance of pause site from that end of the fork that has low BrdU.

    Sample usage:
        < $pause_file python redistribute_sgm_pause_to_inflection_pdf_among_non_paused_forks.py --gaussian-mean $mean\
          --gaussian-sd $sd --fasta $fasta_file
        # in this usage, we set the pause to inflection point pdf to be a gaussian with mean $mean and sd $sd
        # gaussian_mean and gaussian_sd are numbers in bp, pause_file is a pause file,
        # fasta_file is the reference genome. 
        < $pause_file python redistribute_sgm_pause_to_inflection_pdf_among_non_paused_forks.py --gaussian-mean $mean\
          --gaussian-sd $sd --fasta $fasta_file\
          --pause-to-inflection-list $path_to_file_with_pause_to_inflection_point_list\
          --n-bins 20\
          --check-mean-sd-pdf-against-gaussian-with-tol 100
        # in this usage, we set the pause to inflection point pdf using the pause to inflection list in the list file
        # and the number of bins in the options. We also check if the mean and sd of the pause to inflection point pdf
        # matches the gaussian parameters within a tolerance of 100 bp. Other options are as above.
    """

    # get options and data
    df, args = receive_pandas_and_params_from_user([("gaussian-mean", float, 0, "defaults to 0. mean (in bp) if pdf "
                                                                                "is gaussian"),
                                                    ("gaussian-sd", float, 0, "defaults to 0. sd (in bp) if pdf is "
                                                                              "gaussian"),
                                                    ("fasta", str, None, "fasta file of reference genome"),
                                                    ("pause-to-inflection-list", str, "",
                                                     "defaults to empty string. "
                                                     "path to one-column text file with distances between pause sites"
                                                     "and inflection points with the column name rel_pause_loc"),
                                                    ("n-bins", int, 0, "defaults to 0. number of bins in the histogram"
                                                                       "used to construct the pdf from the pause to "
                                                                       "inflection list file. Bins are of uniform width"
                                                                       "and range from min to max of the list."),
                                                    ("check-mean-sd-pdf-against-gaussian-with-tol", float, 0,
                                                     "defaults to 0. If positive, then "
                                                     "check if mean and sd of pause to inflection pdf matches "
                                                     "gaussian parameters within tolerance (bp)")], [],
                                                   desc, sys.stdin,
                                                   required_columns=['detectIndex', 'paramStrNoPause'],
                                                   additional_columns_allowed=True,
                                                   force_tab_separator=True)

    # assign flags based on inputs and perform checks
    is_gaussian_pdf = not os.path.isfile(args.pause_to_inflection_list)
    is_compare_mean_sd = args.check_mean_sd_pdf_against_gaussian_with_tol > 0

    if is_gaussian_pdf and args.gaussian_sd == 0:
        raise ValueError("Please specify a non-zero value for gaussian_sd if you want to use a gaussian pdf")

    if (not is_gaussian_pdf) and args.n_bins < 5:
        raise ValueError("Please specify at least 5 for n_bins if you want to use a pdf from a histogram")

    # set random variable
    if is_gaussian_pdf:
        rv = st.norm(loc=args.gaussian_mean, scale=args.gaussian_sd)
    else:
        data = pd.read_csv(args.pause_to_inflection_list, sep='\s+', comment="#", index_col=None)
        rv = st.rv_histogram(np.histogram(data['rel_pause_loc'], bins=args.n_bins))

    # check if mean and sd of pause to inflection pdf matches gaussian parameters
    if is_compare_mean_sd:
        if abs(rv.mean() - args.gaussian_mean) > args.check_mean_sd_pdf_against_gaussian_with_tol:
            raise ValueError("mean of pause to inflection pdf does not match gaussian mean")
        if abs(rv.std() - args.gaussian_sd) > args.check_mean_sd_pdf_against_gaussian_with_tol:
            raise ValueError("sd of pause to inflection pdf does not match gaussian sd")

    # reject rows where all columns that start with the word "keep" are True i.e. reject valid pauses
    # because we want to redistribute pauses using p.d.f. among non-paused forks
    if df.columns.str.startswith("keep").any():
        df = df[~df.filter(regex='^keep').all(axis=1)]

    # drop all columns except detectIndex, paramStrNoPause
    df = df[['detectIndex', 'paramStrNoPause']]

    # process the detectIndex column using process_fork_index
    df[['__read_id', '__contig', '__start', '__end', '__orn', '__dirn', '__startFork', '__endFork']] \
        = df['detectIndex'].apply(process_fork_index).tolist()

    # set sign of distance depending on fork direction
    df['distSign'] = df['__dirn'].map({'L': -1, 'R': 1})

    # create a column with the inflection point of the noPause sigmoid
    df[['highNoPauseSgm', 'lowNoPauseSgm', 'inflectionPointNoPauseSgm', 'widthNoPauseSgm']] \
        = df['paramStrNoPause'].apply(process_param_str).tolist()

    # print some comments
    print("# pause-inflection point distance pdf was constructed")
    print("# mean of pdf: {}".format(rv.mean()))
    print("# sd of pdf: {}".format(rv.std()))
    print("# fasta file: {}".format(args.fasta))
    if is_gaussian_pdf:
        print("# pdf forced to be gaussian")
    else:
        print("# pdf constructed from histogram of supplied list with {} bins".format(args.n_bins))
        if is_compare_mean_sd:
            print("# supplied mean and sd of pdf for checking: {}, {}".format(args.gaussian_mean, args.gaussian_sd))
            print("# tolerance for mean and sd of pdf against given values: {}".
                  format(args.check_mean_sd_pdf_against_gaussian_with_tol))
        else:
            print("# not performing checks for mean and sd of pdf")

    # print header
    print("detectIndex\tpauseSite\tprobability\tearlyPieceDurn\tlatePieceDurn")

    # loop through each row and calculate probability
    for row_index, row in df.dropna().iterrows():

        # create a list of candidate pause sites, early piece duration, late piece duration
        candidate_pause_sites = [k[1] for k in DNASequenceIterator(row['__contig'], row['__startFork'],
                                                                   row['__endFork'],
                                                                   '+' if row['__orn'] == 'fwd' else '-',
                                                                   args.fasta, 'T')]

        early_piece_durn = [k - row['__startFork'] for k in candidate_pause_sites]
        late_piece_durn = [row['__endFork'] - k for k in candidate_pause_sites]

        # if the fork is left-moving, then switch the duration lists
        if row['__dirn'] == 'L':
            early_piece_durn, late_piece_durn = late_piece_durn, early_piece_durn

        # calculate probability
        probability = rv.pdf([(k - row['inflectionPointNoPauseSgm']) * row['distSign'] for k in candidate_pause_sites])

        # print output
        for k in zip(candidate_pause_sites, probability, early_piece_durn, late_piece_durn):
            print("\t".join([str(entry) for entry in [row['detectIndex'], *k]]))
