import argparse
from Bio import SeqIO, pairwise2

# Set up argument parsing
parser = argparse.ArgumentParser(description="Pairwise comparison of sequences in a FASTA file.")
parser.add_argument("input_fasta", help="Path to the input FASTA file.")
parser.add_argument("output_file", help="Path to the output results file.")
args = parser.parse_args()

# Read all sequences from the input FASTA file
sequences = list(SeqIO.parse(args.input_fasta, "fasta"))

# Filter sequences with read counts >= 2 based on the sequence header (N1-N2-N3-N4 format)
filtered_sequences = []
for seq in sequences:
    # Extract N2 from the header (assuming header format is N1-N2-N3-N4)
    header_parts = seq.id.split('-')
    if len(header_parts) > 1:
        try:
            read_count = int(header_parts[1])
            if read_count >= 2:
                filtered_sequences.append(seq)
        except ValueError:
            # Handle case where N2 is not an integer
            continue

# Function to calculate percent identity
def calculate_percent_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    matches = best_alignment[2]
    aligned_length = max(len(seq1), len(seq2))
    percent_identity = (matches / aligned_length) * 100
    return percent_identity

# Open the output file to write the results
with open(args.output_file, "w") as output_file:
    # Iterate over the filtered sequences
    for seq1 in filtered_sequences:
        # Compare each sequence with all other sequences in the file
        for seq2 in sequences:
            percent_identity = calculate_percent_identity(seq1.seq, seq2.seq)
            # Write the result: header1, header2, percent_identity
            output_file.write(f"{seq1.id}\t{seq2.id}\t{percent_identity:.2f}\n")

print(f"Pairwise comparison completed. Results are saved in '{args.output_file}'.")

