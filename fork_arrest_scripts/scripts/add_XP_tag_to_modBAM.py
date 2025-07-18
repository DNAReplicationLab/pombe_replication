import sys
from DNAscentTools.modBAM_tools_additional import modBAM_add_XP_tag, HeaderLines


if __name__ == "__main__":

    desc = """    Add an XP tag to the end of every .mod.bam line.

    Input: Pipe in modBAM file contents in plain text including headers.

    Output is in plain text and is to stdout. Convert back to .bam if need be using samtools.

    Sample usage:
        samtools view -h sample.mod.bam | python add_XP_tag_to_modBAM.py
    """

    # ensure data is piped in
    if sys.stdin.isatty():
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # object to store header lines
    header = HeaderLines("add_XP_tag_to_modBAM.py", "1.0.0")

    for line in sys.stdin:
        if line.startswith("@"):
            header.add_line(line.strip())
            continue
        else:
            # output header if not already done
            if not header.has_print_action_completed():
                # if program has been run before, then exit
                if header.get_number_of_instances_of_current_program() > 0:
                    raise NotImplementedError("Program has been run before. Exiting.")
                header.add_current_program_to_header()
                header.print_lines()

            # add an XP tag and print modBAM line
            print(modBAM_add_XP_tag(line), end="")
