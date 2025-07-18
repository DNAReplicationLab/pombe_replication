import sys

if __name__ == '__main__':

    # ensure data is piped in
    if sys.stdin.isatty():
        print("Usage: cat sample.detect | python block_mismatched_detect_lines.py")
        print("Goal: When we convert modBAM files with reference-free modification calling and alignment information\n"
              "      to a detect file, we want to remove mismatches i.e. lines where there's a thymidine on the read\n"
              "      but there is some other base on the reference genome. The script should have no effect on a\n"
              "      detect file produced by DNAscent detect or any other reference-based modification calling tool.")
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    rev_flag = False

    for line in sys.stdin:

        stripped_line = line.strip()

        # if a header line, print it
        if stripped_line.startswith("#"):
            print(line, end="")
            continue

        # if a line showing start of a read, mark if it's a reverse read, and print it
        if stripped_line.startswith(">"):
            rev_flag = stripped_line.endswith("rev")
            print(line, end="")
            continue

        # if a line showing a mismatch, skip it, otherwise print it
        first_base_of_kmer = stripped_line[-6]
        last_base_of_kmer = stripped_line[-1]
        if (rev_flag and last_base_of_kmer == "A") or ((not rev_flag) and first_base_of_kmer == "T"):
            print(line, end="")
