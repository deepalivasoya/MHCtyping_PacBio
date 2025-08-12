<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MHC Typing Pipeline</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji";
            line-height: 1.6;
            color: #24292e;
            max-width: 960px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f6f8fa;
        }
        h1, h2, h3 {
            border-bottom: 1px solid #eaecef;
            padding-bottom: 0.3em;
        }
        h1 { font-size: 2em; }
        h2 { font-size: 1.5em; }
        h3 { font-size: 1.25em; }
        pre, code {
            background-color: #f3f4f6;
            border-radius: 6px;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, Courier, monospace;
        }
        pre {
            padding: 16px;
            overflow: auto;
        }
        code {
            padding: 0.2em 0.4em;
        }
        ul {
            padding-left: 2em;
        }
        a {
            color: #0366d6;
            text-decoration: none;
        }
        a:hover {
            text-decoration: underline;
        }
    </style>
</head>
<body>
    <h1>MHC Typing Pipeline for PacBio Data</h1>
    <p>This repository contains a Snakemake pipeline for processing and analyzing PacBio sequencing data of Major Histocompatibility Complex (MHC) genes. The pipeline performs primer trimming, sequence filtering, clustering, and comparison against a reference database to identify potential MHC alleles.</p>

    <h2>Prerequisites</h2>
    <p>Before running the pipeline, you need to have the following software installed:</p>
    <ul>
        <li><strong>Snakemake</strong>: A workflow management system.</li>
        <li><strong>Python 3</strong>: With the following libraries:
            <ul>
                <li><code>Biopython</code></li>
                <li><code>pandas</code></li>
                <li><code>numpy</code></li>
            </ul>
        </li>
        <li><strong>Cutadapt</strong>: A tool for trimming adapters and primers.</li>
        <li><strong>BLAST+</strong>: NCBI's Basic Local Alignment Search Tool.</li>
        <li><strong>Clustal Omega (clustalo)</strong>: A multiple sequence alignment program.</li>
    </ul>
    <p>Or you can use environment.yml file to install conda environment and activate it before running analysis.</p>

    <h2>Setup and Configuration</h2>
    <h3>Directory Structure</h3>
    <p>Organize your project directory as follows:</p>
    <pre><code>.
├── config_&lt;ref&gt;.yaml
├── fasta/
│   ├── &lt;ref&gt;.MHCI.fa
│   └── fullCDS.fa
├── scripts/
│   ├── cluster.py
│   ├── consensus.py
│   ├── fetch_orf.py
│   ├── filter_fasta_length.py
│   ├── flip_reads.py
│   ├── pairwise_comp_clusters.py
│   └── seperate_clusters.py
├── Snakefile
└── fastq/
    └── &lt;Your raw FASTA files here&gt;
</code></pre>
    <ul>
        <li><strong><code>fasta/</code></strong>: This directory should contain reference databases for your species (e.g., <code>Bovine.MHCI.fasta</code>) and a full coding sequence database (<code>fullCDS.fa</code>) used for ORF extraction.</li>
        <li><strong><code>raw_data/</code></strong>: Place your raw PacBio FASTA files here.</li>
    </ul>

    <h3>Configuration File</h3>
    <p>The pipeline's behavior is controlled by a YAML configuration file. Several example files (<code>config_cattle.yaml</code>, <code>config_equine.yaml</code>, <code>config_sheep.yaml</code>) are provided.</p>
    <p>A typical configuration file includes the following parameters:</p>
    <ul>
        <li><code>gene_primer</code>: Forward primer sequence.</li>
        <li><code>gene_primer_rev_comp</code>: Reverse complement of the gene primer.</li>
        <li><code>smart_primer</code>: SMART primer sequence.</li>
        <li><code>umi_seq</code>: Unique molecular identifier sequence.</li>
        <li><code>length_sd</code>: A numerical threshold for filtering sequences based on their length.</li>
        <li><code>miseq_database</code>: Path to the species-specific MHC database.</li>
        <li><code>samples</code>: A list of sample identifiers.</li>
        <li><code>reads</code>: A dictionary mapping sample IDs to their corresponding input FASTA files (required for <code>config_cattle.yaml</code>).</li>
    </ul>
    <p>You must modify the <code>Snakefile</code> to specify the correct configuration file for your analysis.</p>

    <h2>Pipeline Workflow</h2>
    <p>The pipeline is managed by a <code>Snakefile</code> and consists of the following main steps:</p>

    <h3>1. Primer and Adapter Trimming</h3>
    <ul>
        <li><strong><code>trim_gene_primer</code></strong>: Uses Cutadapt to trim gene-specific primers and their reverse complements.</li>
        <li><strong><code>flip_reads</code></strong>: Identifies and reverse-complements reads that were not trimmed in the correct orientation.</li>
        <li><strong><code>trim_smart_primer</code></strong>: Trims the SMART primer sequence.</li>
        <li><strong><code>trim_umi</code></strong>: Trims the UMI sequence from the reads.</li>
    </ul>

    <h3>2. ORF Extraction and Filtering</h3>
    <ul>
        <li><strong><code>get_full_CDS</code></strong>: Uses BLAST to align reads against a <code>fullCDS</code> database and extracts the full open reading frame (ORF).</li>
        <li><strong><code>filter_orf</code></strong>: Filters sequences based on their length, keeping only those within a certain range relative to the mean length.</li>
    </ul>

    <h3>3. Clustering and Comparison</h3>
    <ul>
        <li><strong><code>cluster</code></strong>: Groups sequences into clusters.</li>
        <li><strong><code>mapping</code></strong>: Aligns the clustered sequences against a species-specific MHC database using BLAST.</li>
        <li><strong><code>pairwise</code></strong>: Performs a pairwise comparison of clusters to calculate sequence identity.</li>
        <li><strong><code>cluster_comparison</code></strong>: Separates clusters into "parents" and "children" and generates consensus sequences. This rule also creates a shell script (<code>consensus.sh</code>) to perform multiple sequence alignments and BLAST comparisons of child/singleton consensus sequences against parent clusters.</li>
    </ul>

    <h2>Usage</h2>
    <ol>
        <li><strong>Set the configuration file</strong>: In the <code>Snakefile</code>, change the <code>configfile</code> variable to point to your desired YAML file.
            <pre><code>configfile: "config_sheep.yaml"</code></pre>
        </li>
        <li><strong>Run the pipeline</strong>: Navigate to the top-level directory and execute the following command:
            <pre><code>snakemake --cores &lt;number_of_cores&gt;</code></pre>
            <p>The pipeline will generate all output files in a <code>results/</code> directory.</p>
        </li>
    </ol>

    <h3>Example: Running with Cattle Data</h3>
    <p>If your data corresponds to the <code>config_cattle.yaml</code> file, you would ensure that your <code>Snakefile</code> is configured correctly and your input files are in the <code>raw_data/</code> directory. The pipeline will process the samples listed in the <code>samples</code> and <code>reads</code> sections of the <code>config_cattle.yaml</code> file.</p>
</body>
</html>
