import sys
import argparse

def parse_fai(fai_file):
    """Parse the FASTA index (.fai) file to get contig sizes."""
    contig_sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split the line using whitespace
            parts = line.strip().split()

            # Ensure we have at least two columns for contig name and length
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line in FAI file: {line.strip()}")
                continue

            # Attempt to convert the second part (length) to an integer
            try:
                contig_name = parts[0]  # contig name
                contig_length = int(parts[1])  # contig size
                contig_sizes[contig_name] = contig_length
            except ValueError:
                print(f"Error: Invalid contig size on line: {line.strip()}")
                continue

    return contig_sizes

def parse_origins(origins_file):
    """Parse the origins file, skipping the header row."""
    origins = {}
    with open(origins_file, 'r') as f:
        header = True
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespace

            # Skip empty lines
            if not line:
                continue

            # Process the header
            if header:
                header = False
                continue

            # Split line into parts
            parts = line.split("\t")

            # Ensure we have exactly two parts: contig and position
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line in origins file: {line}")
                continue

            contig = parts[0]
            try:
                position = int(parts[1])
            except ValueError:
                print(f"Warning: Invalid position value: {parts[1]} on line: {line}")
                continue
            
            if contig not in origins:
                origins[contig] = []
            origins[contig].append(position)
    return origins

def count_origins_in_window(contig, start, end, origin_positions):
    """Count origins in the specified window."""
    count = 0
    for origin in origin_positions:
        if start <= origin < end:
            count += 1
    return count

def sliding_window(contig_sizes, origins, window_size, slide_size, output_file):
    """Perform the sliding window analysis and write output to file."""
    with open(output_file, 'w') as out:
        out.write("contig\tstart_window\tend_window\torigin_count\n")
        for contig, size in contig_sizes.items():
            print(f"Processing contig: {contig} of size {size}")  # Debug

            if contig in origins:
                origin_positions = origins[contig]
                print(f"Found {len(origin_positions)} origins for contig {contig}")  # Debug
            else:
                origin_positions = []
                print(f"No origins found for contig {contig}")  # Debug

            for start in range(0, size, slide_size):
                end = min(start + window_size, size)
                origin_count = count_origins_in_window(contig, start, end, origin_positions)

                # Debugging output for the current window
                print(f"Window: {start}-{end}, Origin count: {origin_count}")

                # Write the output for each window
                out.write(f"{contig}\t{start}\t{end}\t{origin_count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count origins per window from genome and origin data.")
    parser.add_argument("fai_file", help="Input FASTA index (.fai) file with contig sizes")
    parser.add_argument("origins_file", help="Input file with origins (tab-separated: contig, position)")
    parser.add_argument("window_size", type=int, help="Window size to use for counting")
    parser.add_argument("slide_size", type=int, help="Sliding window step size")
    parser.add_argument("output_file", help="Output file for writing the results")

    args = parser.parse_args()

    # Parse input files
    contig_sizes = parse_fai(args.fai_file)
    origins = parse_origins(args.origins_file)

    # Perform the sliding window analysis
    sliding_window(contig_sizes, origins, args.window_size, args.slide_size, args.output_file)
