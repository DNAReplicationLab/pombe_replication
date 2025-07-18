import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':

    if len(sys.argv) < 6:
        print('Usage: python plot_pause_count_cd_ho_ranked.py <pause_summary.tsv> <non-zero-label> '
              '<zero-label> <x-axis-label> <output-file-png> <sig-threshold (optional)> <sig-threshold-2 (optional)> '
              ' <publication_flag>')
        print('pause_summary.tsv: a tab-separated file containing our pause summary data')
        print('non-zero-label: the feature label to plot for non-zero pause counts')
        print('zero-label: the feature label to plot for zero pause counts')
        print('NOTE: for the labels above, you can use wildcards e.g. \'Rnh1WTCRAC_Asyn_R1_T.S_no_zeroes\' '
              'matches both TES and TSS as . is a wildcard and stands for any character. If you use wildcards, '
              'then please enclose the label in single quotes on the command line')
        print('x-axis-label: the label for the y-axis of the plot')
        print('output-file-png: the name of the output file to save the plot to')
        print('sig-threshold: optional argument to set the threshold for significance (default is 3), can leave blank')
        print('sig-threshold-lower: optional argument to set another significance mark (default unused, leave blank)')
        print('publication_flag: optional argument to set to True to produce publication formatting, default False')
        sys.exit(1)

    # receive input arguments
    input_file = sys.argv[1]
    non_zero_label = sys.argv[2]
    zero_label = sys.argv[3]
    x_axis_label = sys.argv[4]
    output_file = sys.argv[5]

    if len(sys.argv) > 6 and sys.argv[6] != '':
        sig_threshold = float(sys.argv[6])
    else:
        sig_threshold = 3

    if len(sys.argv) > 7 and sys.argv[7] != '':
        sig_threshold_lower = float(sys.argv[7])
        high_sig = '*'
        low_sig = '†'
    else:
        sig_threshold_lower = None
        high_sig = '*'
        low_sig = ''

    # complain if sig_threshold_lower is set but sig_threshold is not or if sig_threshold_lower is
    # greater than sig_threshold
    if sig_threshold_lower is not None and sig_threshold is None:
        print('Error: sig-threshold-lower is set but sig-threshold is not')
        sys.exit(1)
    if sig_threshold_lower is not None and sig_threshold_lower >= sig_threshold:
        print('Error: sig-threshold-lower is greater than or equal to sig-threshold')
        sys.exit(1)

    is_publication_quality = (len(sys.argv) > 8 and sys.argv[8] == 'True')

    # check that the input file exists and load it
    try:
        df = pd.read_csv(input_file, sep='\s+')
    except FileNotFoundError:
        print('Error: file not found')
        sys.exit(1)

    # check that the output directory exists or make it if it doesn't
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # only retain rows where feature contains the zero label or the non-zero label
    df = df[df['feature'].str.contains(non_zero_label) | df['feature'].str.contains(zero_label)]

    # drop any rows where division is NaN and the feature is not the zero label
    df = df[~(df['division'].isna() & ~df['feature'].str.contains(zero_label))]

    # rename entries in division column by dropping the 'VS_' prefix
    df['division'] = df['division'].str.replace('VS_', 'Q')

    # replace zero or NA values in division column with 'no signal'
    df['division'] = df['division'].str.replace('0', 'no signal')
    df['division'] = df['division'].fillna('no signal')

    # create ratio of observed to expected_sens
    df['ratio'] = df['observed'] / df['expected_sens']

    # create a significant column which contains yes or no based on whether the ratio is sig_threshold standard
    # deviations away from the expected_sens
    df['significant'] = 0
    df.loc[abs(df['observed'] - df['expected_sens']) >= sig_threshold * df['expected_sens_sd'], 'significant'] = 1

    # create another significant column which contains yes or no based on whether the ratio is sig_threshold_lower
    # standard deviations away from the expected_sens
    df['significant_lower'] = 0
    if sig_threshold_lower is not None:
        df.loc[abs(df['observed'] - df['expected_sens']) >= sig_threshold_lower * df['expected_sens_sd'],
               'significant_lower'] = 1

    # drop all columns other than division, feature, relative_direction, ratio, significant, observed, expected_sens,
    # expected_sens_sd, and significant_lower
    df = df[['division', 'feature', 'relative_direction', 'ratio', 'significant', 'observed', 'expected_sens',
             'expected_sens_sd', 'significant_lower']]

    # print the data frame
    print(df)

    # plot bar charts, grouping by relative_direction
    df_list = [df[df['feature'].str.contains('TSS')].groupby(['division', 'relative_direction']),
               df[df['feature'].str.contains('TES')].groupby(['division', 'relative_direction'])]
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    df_list[0]['ratio'].mean().\
        unstack().plot(kind='bar', ax=ax[0], ylim=[0.7, 1.45], yticks=[0.8, 1.0, 1.2, 1.4],
                       color=['#1f77b4', '#ff7f0e'], rot=45, xlabel=x_axis_label,
                       ylabel='Fold enrichment of pauses', hatch='//')
    df_list[1]['ratio'].mean().\
        unstack().plot(kind='bar', ax=ax[1], ylim=[0.7, 1.45], yticks=[0.8, 1.0, 1.2, 1.4],
                       color=['#1f77b4', '#ff7f0e'], rot=45, xlabel=x_axis_label,
                       ylabel='Fold enrichment of pauses', hatch='.')

    # set some plot properties
    for k in range(0, 2):

        # mark significant values
        significant_list = df_list[k]['significant'].mean().unstack().values.flatten('F')
        significant_list_lower = df_list[k]['significant_lower'].mean().unstack().values.flatten('F')
        # NOTE: the flatten option 'F' (column-major) is important above because the default flattens the other way
        # (row-major) and the significant_list will not match up with the patch list below

        for i, p in enumerate(ax[k].patches):
            if significant_list[i] == 1:
                print(p)  # print the patch object for observation/debugging
                ax[k].text(p.get_x() + p.get_width() / 2., p.get_height() + 0.01, high_sig, ha='center', va='bottom',
                           fontsize=15)
            elif significant_list_lower[i] == 1:
                print(p)  # print the patch object for observation/debugging
                ax[k].text(p.get_x() + p.get_width() / 2., p.get_height() + 0.01, low_sig, ha='center', va='bottom',
                           fontsize=15)

        ax[k].spines[['right', 'top']].set_visible(False)
        ax[k].hlines(y=1.0, linewidth=1.5, linestyles='-', xmin=-1, xmax=5, color='#888888')

        if not is_publication_quality:
            ax[k].xaxis.label.set_fontsize(20)
            ax[k].yaxis.label.set_fontsize(20)
            ax[k].legend(title='Relative direction', fontsize=11, title_fontsize=11)
        else:
            # remove axis labels and legend
            ax[k].set_xlabel('')
            ax[k].set_ylabel('')
            ax[k].get_legend().remove()

        ax[k].tick_params(axis='both', which='major', labelsize=15)

    fig.savefig(output_file, dpi=300, facecolor='w', edgecolor='w', bbox_inches='tight', pad_inches=0.1)
