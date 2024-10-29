import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd

def process_fasta(input_fasta, output_fasta, summary_tsv, length_distribution_tsv, threshold):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    seq_lengths = np.array([len(record.seq) for record in sequences])

    #calculate mean and standard deviation
    mean_length = np.mean(seq_lengths)
    sd_length = np.std(seq_lengths)
    min_length = mean_length - threshold
    max_length = mean_length + threshold
    
    selected_sequences = [record for record in sequences if min_length <= len(record.seq) <= max_length]
    removed_sequences = [record for record in sequences if not (min_length <= len(record.seq) <= max_length)]

    with open(output_fasta, "w") as out_fasta:
        for record in selected_sequences:
            out_fasta.write(f">{record.id}\n{str(record.seq)}\n")

    summary_data = {
        "Total sequences": [len(sequences)],
        "Mean length": [mean_length],
        "Standard deviation": [sd_length],
        "Selected sequences": [len(selected_sequences)],
        "Removed sequences": [len(removed_sequences)],
        "Mean length": [mean_length],
        "Threshold": [threshold],
        "Min length considered": [min_length],
        "Max length considered": [max_length]
    }
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    length_distribution = pd.Series(seq_lengths).value_counts().sort_index().reset_index()
    length_distribution.columns = ["Sequence Length", "Frequency"]
    length_distribution.to_csv(length_distribution_tsv, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences based on the mean length and a threshold.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_fasta", help="Output filtered FASTA file")
    parser.add_argument("summary_tsv", help="Output summary TSV file")
    parser.add_argument("length_distribution_tsv", help="Output length distribution TSV file")
    parser.add_argument("--threshold", type=float, required=True, help="Threshold for filtering sequences (mean ± threshold)")
    args = parser.parse_args()

    process_fasta(args.input_fasta, args.output_fasta, args.summary_tsv, args.length_distribution_tsv, args.threshold)

if __name__ == "__main__":
    main()
