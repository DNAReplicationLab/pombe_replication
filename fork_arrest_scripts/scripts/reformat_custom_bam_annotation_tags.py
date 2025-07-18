import re
import argparse
import sys
import uuid


def convert_custom_bam_annotation_tag_to_lists(input_string: str, is_score_float: bool = False) -> dict:
    """
    Converts information in a custom BAM annotation tag to a dictionary of three lists: BED6+1, BED3+1, or plain_text.
    A few examples are shown below.
    The format of the custom tag we use is XT:Z:str_1;str_2;...;str_N where ... means there are 0 or more strings.
    str_i can be of the following formats:
    - str_i = name1,name2,...,nameN where ... means there are 0 or more names
    - str_i = label:contig:start-end,name,score,strand,label:contig:start-end,name,score,strand,...
          where ... means similar entries. The meanings of contig, start, end, name, score, strand are
          the same as in the BED format. The label is a short string that is used to identify the bed file
          from which the entry was obtained (do not use the full filename here but some short stand-in for the file
          name).

    Args:
        input_string: custom BAM annotation tag
        is_score_float: whether the score is a float or not. default is False.

    Returns:
        dictionary that can have none to all of the keys "BED6+1", "BED3+1", or "plain_text", each with a list of lists.

    Examples:
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:gene:chr1:100-200,nom,100,+')
        {'BED6+1': [['chr1', 100, 200, 'nom', 100, '+', 'gene']]}
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:gene:chr1:100-200,,,')
        {'BED3+1': [['chr1', 100, 200, 'gene']]}
        >>> convert_custom_bam_annotation_tag_to_lists('YT:Z')
        {}
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:aa,ba')
        {'plain_text': [['aa'], ['ba']]}
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:dan,ba;;ca,aa;')
        {'plain_text': [['aa'], ['ba'], ['ca'], ['dan']]}
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:gene:chr1:100-200,nom,100,+,gone;'
        ...                                                      'gene2:chr2:400-900,lom,700,-', True)
        {'BED6+1': [['chr1', 100, 200, 'nom', 100.0, '+', 'gene'], ['chr2', 400, 900, 'lom', 700.0, '-', 'gene2']], 'plain_text': [['gone']]}
        >>> convert_custom_bam_annotation_tag_to_lists('XT:Z:gene:chr1:1a0-200,nom,100,+,gone;;;')
        {'plain_text': [['+'], ['100'], ['gene:chr1:1a0-200'], ['gone'], ['nom']]}

    """
    # split the input string into multiple substrings based on the semicolon
    substrings = input_string.strip().split(';')

    # initialize an empty dict to store the output
    output_dict = {"BED6+1": [], "BED3+1": [], "plain_text": []}

    # check that the substrings[0] starts with 'XT:Z:'
    if not substrings[0].startswith('XT:Z:'):
        return {}

    # remove the 'XT:Z:' from the first substring
    substrings[0] = substrings[0][5:]

    # loop through each substring and first split it by the comma
    for i in range(len(substrings)):
        substrings_split = substrings[i].split(',')

        # loop through each substring split by the comma
        # if the entry is not of a suitable format,
        # then add to the output_strings list in the format '# annotate: {entry}'
        # if the entry is of the format, then perform further processing
        # if failure at any point, then add to the output_strings list in the format '# annotate: {entry}'

        j = 0

        while j < len(substrings_split):
            if len(substrings_split[j]) == 0:
                j += 1
                continue
            try:
                label, contig, start_end = substrings_split[j].split(':')
                start, end = start_end.split('-')
                start = int(start)
                end = int(end)
                if not j + 3 < len(substrings_split):
                    raise ValueError
                if re.match(r'^[\x20-\x7e]{1,255}$', substrings_split[j + 1]) and \
                        substrings_split[j + 3] in ['+', '-', '.']:
                    if is_score_float:
                        score = float(substrings_split[j + 2])
                    else:
                        score = int(substrings_split[j + 2])
                    output_dict["BED6+1"].append(
                        [contig, start, end, substrings_split[j + 1], score, substrings_split[j + 3],
                         label])
                    j += 4
                elif substrings_split[j + 1] == substrings_split[j + 2] == substrings_split[j + 3] == '':
                    output_dict["BED3+1"].append([contig, start, end, label])
                    j += 4
                else:
                    raise ValueError
            except ValueError:
                output_dict["plain_text"].append([substrings_split[j]])
                j += 1

    # return the output list with empty lists removed
    return {k: sorted(v) for k, v in output_dict.items() if v}


def process_bam_line_and_output_custom_annotation_tags_as_comments(bam_line: str) -> str:
    """
    Processes a BAM line and outputs the custom annotation tags as comments from XT tags modulated by the XP tag.

    Args:
        bam_line: a data containing line from a bam file

    Returns:
        a string containing potentially several lines corresponding to the annotation information.

    Examples:
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("ff\\tXT:Z:gene:chr1:100-200,nom,100,+")
        '# annotate BED6+1: chr1\\t100\\t200\\tnom\\t100.0\\t+\\tgene'
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("ff\\tXT:Z:gene:chr1:100-200,nom,100,+\\tXP:Z:plot_only_self")
        '# annotate BED6+1: chr1\\t100\\t200\\tnom\\t100.0\\t+\\tgene'
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("5e1a8d61-1730-4356-922b-e031af3ae1ea\\tXT:Z:gene:chr1:100-200,nom,100,+\\tXP:Z:plot_only_self")
        '# annotate BED6+1: chr1\\t100\\t200\\tnom\\t100.0\\t+\\tgene'
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("5e1a8d61-1730-4356-922b-e031af3ae1ea\\tXT:Z:gene:chr1:100-200,5f1a8d61-1730-4356-922b-e031af3ae1ea,100,+\\tXP:Z:plot_only_self")
        '# annotate_no_plot BED6+1: chr1\\t100\\t200\\t5f1a8d61-1730-4356-922b-e031af3ae1ea\\t100.0\\t+\\tgene'
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("5e1a8d61-1730-4356-922b-e031af3ae1ea\\tXT:Z:gene:chr1:100-200,5e1a8d61-1730-4356-922b-e031af3ae1ea,100,+\\tXP:Z:plot_only_self")
        '# annotate BED6+1: chr1\\t100\\t200\\t5e1a8d61-1730-4356-922b-e031af3ae1ea\\t100.0\\t+\\tgene'
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments('YT:Z')
        ''
        >>> process_bam_line_and_output_custom_annotation_tags_as_comments("something\\tsomething\\tXT:Z:aa,ba")
        '# annotate plain_text: aa\\n# annotate plain_text: ba'
    """

    output_string_list = []
    read_id = ""

    name_list = []
    only_plot_self = False

    for k in bam_line.strip().split("\t"):

        # assign read id if not already assigned from the first column
        if read_id == "":
            read_id = k
            continue

        # extract the custom annotation tags
        if k.startswith("XT:Z"):
            for key, list_val in convert_custom_bam_annotation_tag_to_lists(k, True).items():
                for val in list_val:
                    output_string_list.append(f"# annotate {key}: " + "\t".join([str(j) for j in val]))

                    # record the read id from the bed data if available
                    if key == "BED6+1":
                        try:
                            uuid.UUID(val[3])
                            name_list.append(val[3])
                        except ValueError:
                            name_list.append("N/A")
                    else:
                        name_list.append("N/A")

        # if plot_only_self is set, then prepare to filter out other read ids
        only_plot_self = only_plot_self or (k == "XP:Z:plot_only_self")

    if only_plot_self:
        return "\n".join([k if (name_list[i] == "N/A" or name_list[i] == read_id)
                          else k.replace("annotate", "annotate_no_plot") for i, k in enumerate(output_string_list)])
    else:
        return "\n".join(output_string_list)


if __name__ == "__main__":

    desc = """    Get custom annotation tags from modBAM file contents and render as comments (lines starting with '#').

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is to stdout.

    Sample usage:
        samtools view -h sample.mod.bam | python <programName.py>
        samtools view -h sample.mod.bam | python <programName.py> > sample.txt
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # create the comment output
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            output_str = process_bam_line_and_output_custom_annotation_tags_as_comments(line)
            if len(output_str) > 0:
                print(output_str)
