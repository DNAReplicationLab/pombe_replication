import argparse
import sys
from DNAscentTools.modBAM_tools_additional import modBAM_record_windowed_average

if __name__ == "__main__":

    desc = """    Output windowed analogue density per read id from modBAM file.

    Input: Pipe in modBAM file contents in plain text including or excluding headers.

    Output two columns in TSV format: read_id, mean analogue density per window with no header, comments, column names.
    Windows are over bases where modification information is available.
    This is a quick-and-dirty script, so we do not output window positions, and there's no guarantee on windowing
    direction (i.e. left to right or right to left) or offset (i.e. where the first window starts).
    Our only guarantee is window size.
    If you want these additional features, please use other scripts in the repository e.g. get_mean_brdU_window.py.
    Bases where modification probabilities are missing are ignored, not treated as NAs or as zeroes.
    Modbam must follow "?" notation which means that the modification status of a base where modification
    information is missing is 'unclear' as opposed to 'unmodified'. 

    Sample usage:
        # Couple of sample commands are shown below.
        samtools view sample.mod.bam | python <programName.py> --window 300
            # measure BrdU densities in non-overlapping windows of 300 thymidines each
        samtools view sample.mod.bam | python <programName.py> --window 300 --thres 0.10
            # set a threshold of 0.1 different from the default of 0.5
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--window', type=int, required=False,
                        help='(default 3000) number of thymidines per window',
                        default=3000)
    parser.add_argument('--thres', type=float, required=False,
                        help='(default 0.5) threshold above (below) which T is regarded as modified (unmodified)',
                        default=0.5)
    parser.add_argument('--tag', type=str, required=False,
                        help='(default T) ChEBI code or one letter code of base modification',
                        default='T')
    parser.add_argument('--base', type=str, required=False,
                        help='(default T) One letter code of unmodified base',
                        default='T')
    parser.add_argument('--force-missing', action='store_true', required=False,
                        help=('(default False) By default, or by using a "." in the MM tag in the mod bam file,'
                              'authors mandate that missing modification information is interpreted as '
                              'unmodified. Using a ?, authors mandate that missing modification information is '
                              'interpreted as missing. Unfortunately, this practice is not consistent across '
                              'modBAM files. This option allows you to interpret "." or blank as "?" i.e. force the '
                              'interpretation of missing modification information as missing.'),
                        )

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # output read ids of nascent reads
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            read_id, windowed_values = modBAM_record_windowed_average(line.rstrip(), args.window, args.thres,
                                                                      code=args.tag, base=args.base,
                                                                      force_missing=args.force_missing)
            for k in windowed_values:
                print(f"{read_id}\t{k:.6f}")
