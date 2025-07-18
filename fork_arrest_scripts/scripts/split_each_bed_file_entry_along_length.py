import os
import argparse
import json
from validate_bed_format import bed_6_file_from_stdin_to_sorted_pandas_dataframe
from fork_arrest_tools import BedEntryWindowIterator

if __name__ == '__main__':
    desc = """    Align bed file entries by arrow head or tail or along increasing or decreasing genomic coordinates, 
    then split each entry across many files by dividing them according to a given number of bases.
    Arrow head is end column and tail is start column if strand is + and vice versa.
    
    Example: if the bed entry is 'chr1 100 200 +' and the user specifies 40 bases and align by tail, then the first
    file will have the entry 'chr1 100 140 +' and the second file will have the entry 'chr1 140 180 +', and
    the third file will have the entry 'chr1 180 200 +'. This is done for each entry.

    Input: Pipe in a bed file. Also needed are a prefix for the output files and a number of bases to split by.
    If the bed entries are not all of the same length, then the number of entries per file will decrease with the
    numbering of the files as the smaller entries would have been exhausted faster than the larger entries.

    Output: A bunch of bed files with the prefix and a number appended to the end.
    
    Example: cat input.bed | python split_each_bed_file_entry_along_length.py --prefix output --num-bases 1000
        """

    # get options
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prefix', type=str, required=True, help='Prefix for output files.')
    parser.add_argument('--num-bases', type=int, required=True, help='Number of bases (in bp) to split each entry by.')
    parser.add_argument('--align-by', type=str, required=True,
                        choices=['head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref'],
                        help='Align by arrow head or tail or along increasing or decreasing genomic coordinates, '
                             'and proceed from there to the other end.')

    args = parser.parse_args()

    # ensure number of bases is larger than 0
    if args.num_bases < 1:
        raise ValueError('Number of bases must be larger than 0.')

    # get dataframe from bed file in stdin
    df = bed_6_file_from_stdin_to_sorted_pandas_dataframe()

    # get parent directory of prefix
    parent_dir = os.path.dirname(args.prefix)

    # if parent directory doesn't exist, make it
    if not (parent_dir == "" or os.path.exists(parent_dir)):
        os.makedirs(parent_dir, exist_ok=True)

    # store a list of files
    file_list = []

    # iterate through each entry in the dataframe and prepare to write entries to files according to the number of
    # bases requested
    for label, data in df.iterrows():
        count = 0
        for k in BedEntryWindowIterator(data['start'], data['end'], data['strand'], args.num_bases, args.align_by):
            count += 1

            data_copy = data.copy()
            data_copy['start'] = k[0]
            data_copy['end'] = k[1]

            file_name = f'{args.prefix}_{count}.bed'
            file_list.append({"division": str(count), "bed_file": file_name,
                              "data": '\t'.join([str(x) for x in data_copy.values])})

    # first, create all the files
    for file in file_list:
        with open(file['bed_file'], 'w') as f:
            pass

    # then, write the data to the files
    for file in file_list:
        with open(file['bed_file'], 'a') as f:
            f.write(file['data'] + '\n')
        del file['data']

    # print file list to stdout after converting to json
    print(json.dumps([dict(t) for t in {tuple(d.items()) for d in file_list}], indent=4))
