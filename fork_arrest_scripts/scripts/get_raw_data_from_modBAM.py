import pysam
import sys
from clize import run
from user_input_tools import receive_pandas_and_params_from_user
from DNAscentTools.modBAM_tools_additional import process_fork_index, ModBamRecordProcessor


def wrap_get_read_data_from_modBAM(mod_bam_file: str, read_id: str, contig: str, start: int, end: int,
                                   blank_if_fail: bool = False, read_id_replacement: str = "",
                                   get_max: bool = False, base: str = "T", code: str = "T",
                                   return_fwd_seq_coords: bool = False,
                                   move_parallel_fwd_seq: bool = False) -> str:
    """ 
    Returns tabulated raw data from modBAM given a read and coordinates.
    Three tab-sep columns without header read id, coordinate, probability.
    Coordinate can be reference coordinate or fwd seq coordinate in either ascending or descending order
    depending on the input options.
    If the option get_max is set, only the maximum probability is returned, with an NA in the coordinate column.

    Args:
        mod_bam_file: bam file
        read_id: read id of interest
        contig: contig of interest, can be set to "unmapped" or "ref-independent" to specify that the read is unmapped
                or that we want to output modification calls along the forward sequence (i.e. basecalled sequence)
        start: 0-based start position on reference or forward sequence depending on whether
               contig is set to a reference contig or to "unmapped"/"ref-independent"
        end: 0-based end position on reference or forward sequence depending on whether
             contig is set to a reference contig or to "unmapped"/"ref-independent"
        blank_if_fail: (default F) if read is not found, return empty string instead of failing
        read_id_replacement: (default "") in output, replace read_id with read_id_replacement if given
        get_max: (default F) just report the max value instead of all values
        base: base to be considered for modification, default is "T"
        code: code of modification, default is "T". You can use a number here as well.
        return_fwd_seq_coords: (default F) if true, return the forward strand coordinates of the read, not the
            reference coordinates i.e. the modification data as recorded natively instead of mapped to the reference.
            If contig is set to "unmapped" or "ref-independent", this is set to True as only the forward sequence
            coordinates are available.
        move_parallel_fwd_seq: (default F) if true, report values in a direction parallel to the basecalled sequence.
                                defaults to moving parallel to the top strand of the reference sequence.
                                If contig is set to "unmapped" or "ref-independent", this is set to True as only the
                                forward sequence coordinates are available.

    Returns:
        a string containing three columns separated by tabs and row(s) separated by newlines

    """

    # find relevant records
    is_contig_set = False
    if contig == "unmapped":
        alignments = pysam.view("-e", f"qname==\"{read_id}\"", "--include-flags", "UNMAP", mod_bam_file)
        is_contig_set = False
    elif contig == "ref-independent":
        alignments = pysam.view("-e", f"qname==\"{read_id}\"", mod_bam_file)
        is_contig_set = False
    else:
        alignments = pysam.view("-e", f"qname==\"{read_id}\"", mod_bam_file, f"{contig}:{start}-{end}")
        is_contig_set = True

    # threshold below is unused in this script
    mod_bam_parser = ModBamRecordProcessor(0.01, code, allow_non_na_mode=True, base=base)

    # process the data from the bam file
    split_lines = alignments.splitlines()

    # check that only one read with given name is picked up
    index_with_data = [i for i, k in enumerate(split_lines) if k]
    if not (len(index_with_data) == 1):
        if blank_if_fail:
            return ''
        raise NotImplementedError("Cannot find read id or duplicate read ids exist")
    else:
        mod_bam_parser.process_modbam_line(split_lines[index_with_data[0]])

    # return data
    if mod_bam_parser.has_data():
        if is_contig_set:
            read_data = list(filter(lambda y: start <= y.ref_pos < end,
                                    mod_bam_parser.mod_data_to_table(move_parallel_top_ref_strand=
                                                                     not move_parallel_fwd_seq)))
        else:
            read_data = list(filter(lambda y: start <= y.fwd_seq_pos < end,
                                    mod_bam_parser.mod_data_to_table(move_parallel_top_ref_strand=False)))
        if get_max:
            return ((read_id if read_id_replacement == "" else read_id_replacement) + "\tNA\t" +
                    str(max((k.mod_qual for k in read_data), default='NA')))
        else:
            return "\n".join(
                map(
                    lambda x: "\t".join((read_id if read_id_replacement == "" else read_id_replacement,
                                         str(x.ref_pos) if (is_contig_set and not return_fwd_seq_coords)
                                         else str(x.fwd_seq_pos), str(x.mod_qual))),
                    read_data
                )
            )
    else:
        if blank_if_fail:
            return ''
        raise NotImplementedError("No data found")


def piped_regions(mod_bam_file: str, *, alt_read_id_column: bool = False,
                  infer_contig_read_id: bool = False,
                  fork_info_in_alt_read_id: bool = False,
                  get_max: bool = False,
                  base: str = "T", code: str = "T") -> str:
    """
    Pipe in 4 columns with header read_id contig start end, separated by spaces or tabs.

    Args:
        mod_bam_file: Location of mod bam file
        alt_read_id_column: (default F) if true, an additional 'alt_read_id' column is present,
            and program replaces read_id with alt_read_id in the output
        infer_contig_read_id: (default F) deprecated. Use fork_info_in_alt_read_id instead if you
            want to extract contig and read_id from alt_read. Retained for backward compatibility.
            Option will throw an error if used.
        fork_info_in_alt_read_id: (default F) if true, 'alt_read_id' column must have the format
            "readID_contig_start_end_orn_dirn_startFork_endFork". readID, contig, startFork, endFork
            are extracted and used to get raw data from the modBAM file. Automatically
            sets alt_read_id_column to True.
        get_max: (default F) report only the maximum probability per read_id. If set, the position column is returned
            as all NAs.
        base: (default "T") base to be considered for modification
        code: (default "T") code of modification. You can use a number here as well.

    Returns:
        String of three columns and separated by tabs of read_id, position, probability

    """

    # receive input from user
    df, _ = receive_pandas_and_params_from_user([], [], '', sys.stdin,
                                                required_columns=[])

    if fork_info_in_alt_read_id:
        df[['read_id', 'contig', '__start', '__end', '__orn', '__dirn', 'start', 'end']] \
            = df['alt_read_id'].apply(process_fork_index).tolist()
    elif infer_contig_read_id:
        raise NotImplementedError("infer_contig_read_id is deprecated. Use fork_info_in_alt_read_id instead.")
    elif not alt_read_id_column:
        df['alt_read_id'] = ""

    return "\n".join(k for k in
                     (wrap_get_read_data_from_modBAM(mod_bam_file, row[0], row[1], row[2], row[3], blank_if_fail=True,
                                                     read_id_replacement=row[4], get_max=get_max, base=base, code=code)
                      for row in zip(df['read_id'], df['contig'], df['start'], df['end'], df['alt_read_id']))
                     if len(k) > 0)


if __name__ == '__main__':
    # to get help, run python program_name.py --help
    run(wrap_get_read_data_from_modBAM, alt=piped_regions)
