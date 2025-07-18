import pandas as pd
import sys
from io import StringIO

if __name__ == "__main__":

    # ensure data is piped in
    if sys.stdin.isatty():
        print("Goal: Create a LaTeX document with all plots grouped by feature and plot_file")
        print("Usage: python collect_plots_per_dataset_and_feature_to_latex.py < data.tsv > latex_document.tex")
        print("data.tsv should be a tab-separated file with columns 'folder', 'plot_file', 'dataset', 'feature'")
        print("Any comments in data.tsv should start with '#' and are printed upfront in the LaTeX document")
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # read all lines from stdin
    lines = sys.stdin.readlines()

    # extract the comments from the lines
    comments = [line for line in lines if line.startswith('#')]

    # Read the TSV dataset using pandas
    data = pd.read_csv(StringIO("\n".join(lines)), delimiter='\t', comment='#')

    # check that the required columns are present
    required_columns = ['folder', 'plot_file', 'dataset', 'feature']
    for column in required_columns:
        if column not in data.columns:
            raise ValueError(f"Column {column} not found in data.tsv")

    # check that there is at least one row
    if data.shape[0] == 0:
        raise ValueError("data tsv is empty")

    # Start the LaTeX document
    latex_content = '''
\\documentclass{article}
\\usepackage{graphicx}
\\usepackage{hyperref}
\\usepackage[a4paper, margin=0.5in]{geometry}
\\title{Compare pause enrichment across different datasets}
\\date{\\today}
\\begin{document}
\\maketitle
\\tableofcontents
\\pagebreak
\\section{Comments in input file}
'''

    # Add the comments to the LaTeX document
    for comment in comments:
        latex_content += comment[1:].replace('_', '\\_')  # Remove the leading '#'

    latex_content += '\\pagebreak\n'  # Start a new page after the comments

    # Group by 'feature'
    for feature, feature_group in data.groupby('feature'):
        latex_content += '\\section{' + feature.replace('_', '\\_') + '}\n'

        # Group by 'plot_file' within each 'feature' group
        for plot_file, plot_file_group in feature_group.groupby('plot_file'):
            latex_content += '\\subsection{' + plot_file.replace('_', '\\_') + '}\n'

            # Add subsubsection for each 'dataset' and the associated image
            for index, row in plot_file_group.iterrows():
                latex_content += '\\subsubsection{' + row['dataset'].replace('_', '\\_') + '}\n'
                latex_content += '\\includegraphics[width=\\linewidth]{' + row['folder'] + '/' + row[
                    'plot_file'] + '}\n'

            latex_content += '\\newpage\n'  # Start a new page for each 'plot_file'

    # End the LaTeX document
    latex_content += '\\end{document}'

    # Write the LaTeX content to stdout
    sys.stdout.write(latex_content)
