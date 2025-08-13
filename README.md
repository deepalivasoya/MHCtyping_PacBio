# **MHC Typing Pipeline for PacBio Data**

This repository contains a Snakemake pipeline for processing and analyzing PacBio sequencing data of Major Histocompatibility Complex (MHC) genes. The pipeline performs primer trimming, sequence filtering, clustering, and comparison against a reference database to identify potential MHC alleles.

## **Prerequisites**

Before running the pipeline, you need to have the following software installed:

* **Snakemake**: A workflow management system.  
* **Python 3**: With the following libraries:  
  * Biopython  
  * pandas  
  * numpy  
* **Cutadapt**: A tool for trimming adapters and primers.  
* **BLAST+**: NCBI's Basic Local Alignment Search Tool.  
* **Clustal Omega (clustalo)**: A multiple sequence alignment program.

Alternatively, you can use the environment.yml file to create and activate a Conda environment before running the analysis.

## **Setup and Configuration**

### **Directory Structure**

Organize your project directory as follows:

.  
├── config\_\<ref\>.yaml  
├── fasta/  
│   ├── \<ref\>.MHCI.fa  
│   └── fullCDS.fa  
├── scripts/  
│   ├── cluster.py  
│   ├── consensus.py  
│   ├── fetch\_orf.py  
│   ├── filter\_fasta\_length.py  
│   ├── flip\_reads.py  
│   ├── pairwise\_comp\_clusters.py  
│   └── seperate\_clusters.py  
├── Snakefile  
└── fastq/  
    └── \<Your raw FASTA files here\>

* **fasta/**: This directory should contain reference databases for your species (e.g., Bovine.MHCI.fasta) and a full coding sequence database (fullCDS.fa) used for ORF extraction.  
* **fastq/**: Place your raw PacBio FASTA files here.

### **Configuration File**

The pipeline's behavior is controlled by a YAML configuration file. Several example files (config\_cattle.yaml, config\_equine.yaml, config\_sheep.yaml) are provided.

A typical configuration file includes the following parameters:

* gene\_primer: Forward primer sequence.  
* gene\_primer\_rev\_comp: Reverse complement of the gene primer.  
* smart\_primer: SMART primer sequence.  
* umi\_seq: Unique molecular identifier sequence.  
* length\_sd: A numerical threshold for filtering sequences based on their length.  
* miseq\_database: Path to the species-specific MHC database.  
* samples: A list of sample identifiers.  
* reads: A dictionary mapping sample IDs to their corresponding input FASTA files (required for config\_cattle.yaml).

You must modify the Snakefile to specify the correct configuration file for your analysis.

## **Pipeline Workflow**

The pipeline is managed by a Snakefile and consists of the following main steps:

### **1\. Primer and Adapter Trimming**

* **trim\_gene\_primer**: Uses Cutadapt to trim gene-specific primers and their reverse complements.  
* **flip\_reads**: Identifies and reverse-complements reads that were not trimmed in the correct orientation.  
* **trim\_smart\_primer**: Trims the SMART primer sequence.  
* **trim\_umi**: Trims the UMI sequence from the reads.

### **2\. ORF Extraction and Filtering**

* **get\_full\_CDS**: Uses BLAST to align reads against a fullCDS database and extracts the full open reading frame (ORF).  
* **filter\_orf**: Filters sequences based on their length, keeping only those within a certain range relative to the mean length.

### **3\. Clustering and Comparison**

* **cluster**: Groups sequences into clusters.  
* **mapping**: Aligns the clustered sequences against a species-specific MHC database using BLAST.  
* **pairwise**: Performs a pairwise comparison of clusters to calculate sequence identity.  
* **cluster\_comparison**: Separates clusters into "parents" and "children" and generates consensus sequences. This rule also creates a shell script (consensus.sh) to perform multiple sequence alignments and BLAST comparisons of child/singleton consensus sequences against parent clusters.

## **Usage**

1. **Set the configuration file**: In the Snakefile, change the configfile variable to point to your desired YAML file.  
   configfile: "config\_sheep.yaml"

2. **Run the pipeline**: Navigate to the top-level directory and execute the following command:  
   snakemake \--cores \<number\_of\_cores\>

   The pipeline will generate all output files in a results/ directory.

### **Example: Running with Cattle Data**

If your data corresponds to the config\_cattle.yaml file, you would ensure that your Snakefile is configured correctly and your input files are in the fastq/ directory. The pipeline will process the samples listed in the samples and reads sections of the config\_cattle.yaml file.
