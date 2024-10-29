import sys

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def read_fasta(file_path):
    sequences = {}
    header = None
    sequence = []
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.upper())
        if header:
            sequences[header] = ''.join(sequence)
    
    return sequences

def write_fasta(sequences, output_file):
    #to write sequences to a new FASTA file
    with open(output_file, 'w') as out_file:
        for header, seq in sequences.items():
            out_file.write(f">{header}\n")
            out_file.write(seq + '\n')

def process_fasta(fasta_file, ids_file, output_file):
    #Read list of sequence IDs
    with open(ids_file, 'r') as f:
        ids_to_keep = {line.strip() for line in f if line.strip()}
    
    sequences = read_fasta(fasta_file)
    for header in sequences:
        if header in ids_to_keep:
            sequences[header] = reverse_complement(sequences[header])
    
    write_fasta(sequences, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python reverse_complement_fasta.py <fasta_file> <ids_file> <output_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    ids_file = sys.argv[2]
    output_file = sys.argv[3]
    
    process_fasta(fasta_file, ids_file, output_file)

