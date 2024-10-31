import sys
from Bio import SeqIO
import os

in_pairwise = sys.argv[1]
in_clusters = sys.argv[2]
in_blast = sys.argv[3]
out_prefix = sys.argv[4]
out_folder = sys.argv[5]

out_summary = f"{out_prefix}_cluster_summary.tsv"
out_parent = f"{out_prefix}_parents_clusters.fa"
out_discarded = f"{out_prefix}_discarded_clusters.fa"

parents = {}
cluster = {}
blast = {}
order = {}
reads = {}

# Read the pairwise file
with open(in_pairwise) as f:
    for line in f:
        words = line.strip().split("\t")
        parent = words[0]
        child = words[1]
        if words[2] == '100.00' and child == parent and parent not in cluster:
            parents[parent] = 0
            cluster[parent] = parent
        elif float(words[2]) >= 98 and child not in cluster and parent in parents:
            cluster[child] = parent
            parents[parent] += 1

# Read the BLAST file
with open(in_blast) as f:
    for line in f:
        words = line.strip().split("\t")
        blast[words[0]] = f"{words[1]}|{words[2]}|{words[3]}"

# Write the summary file
total_clusters = 0
total_reads = 0
total_parent_reads = 0
total_singletons = 0
count = 0

with open(f"{out_summary}", "w") as out_summary_file:
    for parent in sorted(parents, key=parents.get, reverse=True):
        count += 1
        info = parent.split("-")
        total_parent_reads += int(info[1])
        order[parent] = count
        count_reads = 0
        singleton_reads = 0
        total_clusters = total_clusters + parents[parent]
        for child in cluster:
            if cluster[child] == parent and child != parent:
                info = child.split("-")
                if info[1] == "1":
                    singleton_reads += 1
                    total_singletons += 1
                else:
                    # print(child)
                    total_reads += int(info[1])
                    count_reads += int(info[1])
        reads[parent] = count_reads
        out_summary_file.write(f"{count}\t{parent}\t{parents[parent]}\t{count_reads}\t{singleton_reads}\t{blast.get(parent, '')}\n")

total_parents = len(parents)
print(f"Total parents: {total_parents}")
print(f"Total parent reads: {total_parent_reads}")
print(f"Total clusters: {total_clusters}")
print(f"Total cluster reads: {total_reads}")
print(f"Total Singletons: {total_singletons}")

os.makedirs(out_folder, exist_ok=True)

# Write the parent and discarded files
with open(f"{out_folder}/{out_parent}", "w") as out_parent_file, open(f"{out_folder}/{out_discarded}", "w") as out_dis_file:
    count_dis = 0
    reads_dis = 0

    # Parse the clusters fasta file using Biopython
    with open(in_clusters) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id = record.id
            seq = str(record.seq)
            info = id.split("-")
            
            if id in parents:
                header = f"{order[id]}|{id}|{blast.get(id, '')}|{parents[id]}|{reads[id]}"
                out_parent_file.write(f">{header}\n{seq}\n")
            elif id in cluster:
                parent = cluster[id]
                count = order[parent]
                if info[1] == "1":
                    with open(f"{out_folder}/{count}_singleton.fa", "a") as singleton_file:
                        singleton_file.write(f">{id}|{blast.get(id, '')}\n{seq}\n")
                else:
                    with open(f"{out_folder}/{count}_child.fa", "a") as child_file:
                        for _ in range(int(info[1])):
                            child_file.write(f">{id}|{blast.get(id, '')}\n{seq}\n")
                    with open(f"{out_folder}/{count}_child_clusters.fa", "a") as child_file:
                            child_file.write(f">{id}|{blast.get(id, '')}\n{seq}\n")
            else:
                out_dis_file.write(f">{id}|{blast.get(id, '')}\n{seq}\n")
                count_dis += 1
                reads_dis += int(info[1])

print(f"Total Discarded: {count_dis}")
print(f"Total Discarded reads: {reads_dis}")

with open(f"{out_folder}/consensus.sh", "w") as consensus:
    consensus.write(f"makeblastdb -dbtype nucl -in {out_folder}/{out_prefix}_parents_clusters.fa\n")
    for cluster in parents:
        count = order[cluster]
        child = f"{out_folder}/{count}_child.fa"
        singleton = f"{out_folder}/{count}_singleton.fa"
        if os.path.exists(child):
            consensus.write(f"clustalo -i {out_folder}/{count}_child.fa -o {out_folder}/{count}_child.aln --outfmt clu --threads 8 --force\n")
            consensus.write(f"python scripts/consensus.py {out_folder}/{count}_child.aln {out_folder}/{count}_child.cons.fa {count}_child_cons\n")
            consensus.write(f"blastn -db {out_folder}/{out_prefix}_parents_clusters.fa -query {out_folder}/{count}_child.cons.fa -outfmt \"6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen\" -out {out_folder}/{count}_child.cons.blast\n")
        if os.path.exists(singleton):
            consensus.write(f"clustalo -i {out_folder}/{count}_singleton.fa -o {out_folder}/{count}_singleton.aln --outfmt clu --threads 8 --force\n")
            consensus.write(f"python scripts/consensus.py {out_folder}/{count}_singleton.aln {out_folder}/{count}_singleton.cons.fa {count}_singleton_cons\n")
            consensus.write(f"blastn -db {out_folder}/{out_prefix}_parents_clusters.fa -query {out_folder}/{count}_singleton.cons.fa -outfmt \"6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen\" -out {out_folder}/{count}_singleton.cons.blast\n")
                    
