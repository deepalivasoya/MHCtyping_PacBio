configfile: "config_cattle.yaml"
# print(config['samples'][0])

rule main:
  input: expand("results/{sample}/clusters/{sample}_parents_clusters.fa", sample=config["samples"])
 
rule trim_gene_primer:
  input: lambda wildcards: "reads/{}.format(config['reads'][wildcards.sample]),
  output: "results/{sample}/cutadapt1.log"
  params:
   prefix = "{sample}",
   gene_primer = config["gene_primer"],
   gene_primer_rev_comp =  config["gene_primer_rev_comp"]
  shell:
    "cutadapt -b {params.gene_primer} -b {params.gene_primer_rev_comp} --report full --info-file results/{params.prefix}/cutadapt1.info --untrimmed-output results/{params.prefix}/cutadapt1.untrimmed -o results/{params.prefix}/cutadapt1.fa {input} > {output}"


rule flip_reads:
  input: "results/{sample}/cutadapt1.log"
  output: "results/{sample}/flipped.fa"
  params:
   prefix = "{sample}",
   cmd = r"""-F $'\t' '{if ($3 == 0) print $1}'"""
  shell:
    r"""
      awk {params.cmd} results/{params.prefix}/cutadapt1.info > results/{params.prefix}/flip_seq_ids.txt 
      python scripts/flip_reads.py results/{params.prefix}/cutadapt1.fa results/{params.prefix}/flip_seq_ids.txt results/{params.prefix}/flipped.fa
    """

rule trim_smart_primer:
  input: "results/{sample}/flipped.fa"
  output: "results/{sample}/cutadapt2.log"
  params:
   prefix = "{sample}",
   smart_primer = config["smart_primer"]
  shell: "cutadapt -g {params.smart_primer} --report full --info-file results/{params.prefix}/cutadapt2.info --untrimmed-output results/{params.prefix}/cutadapt2.untrimmed.fa -o results/{params.prefix}/cutadapt2.fa results/{params.prefix}/flipped.fa > results/{params.prefix}/cutadapt2.log"

rule trim_umi:
  input: "results/{sample}/cutadapt2.log"
  output: "results/{sample}/cutadapt3.log"
  params:
   prefix = "{sample}",
   umi_seq = config["umi_seq"]
  shell: "cutadapt -g {params.umi_seq} --report full --info-file results/{params.prefix}/cutadapt3.info --untrimmed-output results/{params.prefix}/cutadapt3.untrimmed.fa -o results/{params.prefix}/cutadapt3.fa results/{params.prefix}/cutadapt2.fa > results/{params.prefix}/cutadapt3.log"

rule get_full_CDS:
  input: "results/{sample}/cutadapt3.log"
  output: "results/{sample}/orf.fa"
  params:
   prefix = "{sample}",
   cmd = r"""'{if ($9 == 1 || $9 == 2) print $0}'"""
  shell:
    r"""
      blastn -db fasta/fullCDS.fa -query results/{params.prefix}/cutadapt3.fa -outfmt "6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen" -out results/{params.prefix}/fullCDA.blast
      awk {params.cmd} results/{params.prefix}/fullCDA.blast | cut -f1,6 > results/{params.prefix}/orf.location
      python scripts/fetch_orf.py results/{params.prefix}/cutadapt3.fa results/{params.prefix}/orf.location results/{params.prefix}/orf.fa
    """

rule filter_orf:
  input: "results/{sample}/orf.fa"
  output: "results/{sample}/orf.lenSelected.fa"
  params:
   prefix = "{sample}",
   len_threshold = config["length_sd"]
  shell: "python scripts/filter_fasta_length.py {input} {output} results/{params.prefix}/orf_seq_len_report.tsv results/{params.prefix}/orf_length_distribution.tsv --threshold {params.len_threshold}"


rule cluster:
  input: "results/{sample}/orf.lenSelected.fa"
  output: "results/{sample}/clusters.fa"
  shell: 
    "python scripts/cluster.py {input} {output}"

rule mapping:
  input: "results/{sample}/clusters.fa"
  output: "results/{sample}/miseq_all_clusters.blast"
  params:
   prefix = "{sample}",
   database = config["miseq_database"]
  shell:
    r"""
     blastn -db {params.database} -query {input} -outfmt "6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen" -out {output} -max_target_seqs 1
    """

rule pairwise:
  input: "results/{sample}/miseq_all_clusters.blast"
  output: "results/{sample}/clusters.pairwise.comp.tsv"
  params:
   prefix = "{sample}"
  shell:
    r"""
      python scripts/pairwise_comp_clusters.py results/{params.prefix}/miseq.fa results/{params.prefix}/miseq.pairwise_comp.tsv
      python scripts/pairwise_comp_clusters.py results/{params.prefix}/clusters.fa {output}
    """

rule cluster_comparison:
  input: "results/{sample}/miseq_all_clusters.blast"
  output: "results/{sample}/clusters/{sample}_parents_clusters.fa"
  params:
   prefix = "{sample}",
   folder = "results/{sample}/clusters"
  shell:
    r"""
      if [ -d {params.folder} ]; then
        rm -r {params.folder}
      else
        mkdir {params.folder}
      fi
      if [ -f {params.folder}/consensus.sh ]; then
        rm {params.folder}/consensus.sh
      fi
      mkdir {params.folder}
      python scripts/seperate_clusters.py results/{params.prefix}/clusters.pairwise.comp.tsv results/{params.prefix}/clusters.fa results/{params.prefix}/miseq_all_clusters.blast {params.prefix} {params.folder} > results/{params.prefix}/seperate_clusters.log
      makeblastdb -in {params.folder}/{params.prefix}_parents_clusters.fa -dbtype nucl
      bash {params.folder}/consensus.sh
    """
