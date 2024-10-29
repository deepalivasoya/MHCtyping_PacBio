import sys

def parse_fasta(file_path):
    fasta_dict = {}
    header = None
    sequence = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    fasta_dict[header] = ''.join(sequence)
                header = line.split()[0][1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            fasta_dict[header] = ''.join(sequence)
    return fasta_dict

def parse_positions(file_path):
    positions = []
    with open(file_path, 'r') as file:
        for line in file:
            header, start = line.strip().split()
            positions.append((header, int(start)))
    return positions

def fetch_trimmed_sequences(fasta_dict, positions, output_file):
    with open(output_file, 'w') as output:
        for header, start in positions:
            if header in fasta_dict:
                sequence = fasta_dict[header]
                trimmed_sequence = sequence[start-1:]
                if trimmed_sequence:
                    output.write(f">{header}\n")
                    output.write(f"{trimmed_sequence}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <file_A.fasta> <file_B.txt> <output.fasta>")
        sys.exit(1)

    file_A = sys.argv[1]
    file_B = sys.argv[2]
    output_file = sys.argv[3]
    fasta_dict = parse_fasta(file_A)
    positions = parse_positions(file_B)
    fetch_trimmed_sequences(fasta_dict, positions, output_file)

if __name__ == "__main__":
    main()

