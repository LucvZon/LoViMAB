# LoViMAB
Long-Read Viral Metagenome Assembly Benchmark (LoViMAB) is a scalable and reproducible Snakemake workflow designed to benchmark the performance of various long-read assemblers on viral metagenomic Nanopore data.

Assembling complete viral genomes from complex metagenomes is a significant challenge. This pipeline addresses this by implementing a "semi-targeted" approach:

1. **Targeted Read Filtering**: Pre-filters raw Nanopore reads using a custom protein database (e.g., all human viruses) to enrich for sequences of interest.
2. **Comparative Assembly**: Assembles the filtered reads using multiple state-of-the-art assemblers.
3. **In-depth Quality Assessment**: Performs comprehensive, reference-free ([CheckV](https://bitbucket.org/berkeleylab/checkv)) and reference-based ([QUAST](https://github.com/ablab/quast)) quality control to generate detailed comparison reports.

This workflow is designed to help researchers evaluate which assembly tool is most appropriate for their data, especially in scenarios where reliable reference genomes are unavailable or highly diverse.

Current supported assemblers:
- [metaFlye](https://github.com/mikolmogorov/Flye)
- [Canu](https://github.com/marbl/canu)
- [Raven](https://github.com/lbcb-sci/raven)
- [PenguiN](https://github.com/soedinglab/plass)

# How to run

```
# Clone repo
git clone https://github.com/LucvZon/LoViMAB.git

cd LoViMAB

# Install environment
conda env create -f env/environment.yml

# Activate environment
conda activate

snakemake --snakefile workflow/snakefile.smk --cores 24
```

# Setup configuration

LoViMAB relies on **config/lovimab.yml** to find input fastq data, databases, references and more.

Key sections to adjust:

- **Sample Configuration**: Choose a sample name and supply a path to its corresponding fastq.gz files.
- **External Database and Script Paths**: Set the diamond and checkv databases.
- **Virus Selection**: Fill in a species name and its corresponding FASTA reference file.

Optional sections:

- **Assembler Selection**: Turn assemblers on/off.
- **General Parameters**: Adjust number of threads and specific tool parameters

```
# -------------------------------------------------------------------
#          CONFIGURATION FOR TARGETED ASSEMBLY WORKFLOW
# -------------------------------------------------------------------
 
# --- Sample Configuration ---
# List your samples here. The key is the sample name (e.g., "Sample1"),
# and the value is the path and prefix of the input fastq files.
# The workflow will find all *.fastq.gz files in your path.
samples:
  Sample1: "/path/to/input/reads"
  # Sample2: "path/to/input/reads"
 
# --- External Database and Script Paths ---
# Fill in the absolute paths to your databases and scripts.
paths:
  # Database for initial read classification
  diamond_db: "/path/to/diamond/database.dmnd"
  # Database for CheckV
  checkv_db: "/path/to/checkv/"
  # Path to diamond post-processing script
  post_process_script: "/path/to/post_process_diamond_v1.0.py"
 
# --- Assembler Selection ---
# Set to 'true' to run an assembler, or 'false' to skip it.
assemblers:
  metaflye: true
  penguin: true
  raven: true
  canu: true
  # Add other assemblers if needed
 
# --- General Parameters ---
params:
  threads: 24
  # Minimum contig length for PenguiN
  penguin_min_contig_len: 500
  # Minimum sequence identity for PenguiN overlaps
  penguin_min_seq_id: 0.95
  # Set genome size for canu
  canu_genome_size: 200k
  # Set DIAMOND BLASTX sensitivity for read filtering
  # options: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive
  diamond_sensitivity: sensitive

# --- Virus Selection ---
# Add any viruses you want to specifically deep-dive on
viruses_of_interest:
  # Key: A species name for the virus
  # Value: Path to the reference FASTA for that virus
  Monkeypox_virus: "path/to/reference/genome"
  Influenza_A_virus: "path/to/reference/genome"
```
