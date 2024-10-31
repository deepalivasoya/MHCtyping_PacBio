import sys
from collections import defaultdict

def parse_fasta(filename):
    sequences = []
    header = None
    sequence = []
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, ''.join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            sequences.append((header, ''.join(sequence)))
    return sequences

def write_fasta(clusters, output_file):
    total_sequences = sum(cluster['count'] for cluster in clusters)
    
    with open(output_file, 'w') as file:
        for idx, cluster in enumerate(clusters, 1):
            seq_count = cluster['count']
            percent = round((seq_count / total_sequences) * 100, 2)
            length = len(cluster['sequence'])
            header = f"{idx}-{seq_count}-{percent}-{length}"
            file.write(f">{header}\n{cluster['sequence']}\n")

def write_report(clusters, report_file):
    total_clusters = len(clusters)
    single_sequence_clusters = sum(1 for cluster in clusters if cluster['count'] == 1)
    
    with open(report_file, 'w') as file:
        file.write(f"{total_clusters}\t{single_sequence_clusters}\n")

def cluster_sequences(sequences):
    seq_dict = defaultdict(int)
    
    for _, seq in sequences:
        seq_dict[seq] += 1
    
    clusters = [{'sequence': seq, 'count': count} for seq, count in seq_dict.items()]
    clusters.sort(key=lambda x: -x['count'])
    
    return clusters

def main(input_file, output_fasta, report_file):
    sequences = parse_fasta(input_file)
    clusters = cluster_sequences(sequences)
    write_fasta(clusters, output_fasta)
    write_report(clusters, report_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python cluster_fasta.py <input_fasta> <output_fasta> <report_file>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    report_file = sys.argv[3]
    
    main(input_fasta, output_fasta, report_file)

