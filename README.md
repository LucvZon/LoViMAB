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

Clone the repo and install the required conda environment:

```
# Clone repo
git clone https://github.com/LucvZon/LoViMAB.git

cd LoViMAB

# Install environment
conda env create -f env/environment.yml

# Activate environment
conda activate lovimab
```

## Setup configuration

LoViMAB relies on `config/lovimab.yml` to find input fastq data, databases, references and more.

Key sections to adjust:

- `Sample Configuration`: Choose a sample name and supply a path to its corresponding fastq.gz files.
- `External Databases`: Set the paths for the diamond and checkv databases.
- `Virus Selection`: Fill in a species name and its corresponding FASTA reference file.

Optional sections:

- `Assembler Selection`: Turn assemblers on/off.
- `General Parameters`: Adjust number of threads and specific tool parameters

**config/lovimab.yml**:
```
# -------------------------------------------------------------------
#          CONFIGURATION FOR TARGETED ASSEMBLY WORKFLOW
# -------------------------------------------------------------------
 
# --- Sample Configuration ---
# List your samples here. The key is the sample name (e.g., "Sample1"),
# and the value is the path and prefix of the input fastq files.
# The workflow will find all *.fastq.gz files in your path.
samples:
  Sample1: "/path/to/input/reads/"
  # Sample2: "path/to/input/reads/"
 
# --- External Databases ---
# Fill in the absolute paths to your databases and scripts.
paths:
  # Database for initial read classification
  diamond_db: "/path/to/diamond/database.dmnd"
  # Database for CheckV
  checkv_db: "/path/to/checkv/"
 
# --- Assembler Selection ---
# Set to 'true' to run an assembler, or 'false' to skip it.
assemblers:
  metaflye: true
  penguin: true
  raven: true
  canu: true
  myloasm: true
  metamdbg: true
  wtdbg2: true
  shasta: true
  miniasm: true
  # Add other assemblers if needed
 
# --- General Parameters ---
params:
  threads: 24

  # --- PENGUIN ---
  # Minimum contig length
  penguin_min_contig_len: 500
  # Minimum sequence identity for overlaps
  penguin_min_seq_id: 0.95

  # --- CANU ---
  # Set genome size for canu
  canu_genome_size: 200k

  # --- MYLOASM ---
  myloasm_min_reads: 10
  myloasm_min_overlap: 250

  # --- METAMDBG ---
  metamdbg_min_overlap: 250
  metamdbg_min_seq_id: 0.95

  # --- WTDBG2 ---
  wtdbg2_min_read_len: 500
  wtdbg2_min_contig_len: 1000

  # --- SHASTA ---
  shasta_min_read_len: 500

  # --- MINIASM ---
  miniasm_min_overlap: 400

  # --- DIAMOND ---
  # Set DIAMOND BLASTX sensitivity for read filtering
  # options: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive
  diamond_sensitivity: sensitive

# --- Virus Selection ---
# Add any viruses you want to specifically deep-dive on
viruses_of_interest:
  # Key: A species name for the virus
  # Value: Path to the reference FASTA for that virus
  Virus_name_1: "path/to/reference/genome"
  Virus_name_2: "path/to/reference/genome"
```

If your conda environment and config/lovibam.yml are ready, then run the workflow:
```
snakemake --snakefile workflow/snakefile.smk --cores 24
```
Adjust --cores 24 if necessary. 

# Workflow

### 1. Read Preparation & QC
- **Input**: Raw FASTQ files for a sample.
- **Process**:
	- -> Merge all raw reads.
	- -> Trim adapters and low-quality ends ([cutadapt](https://github.com/marcelm/cutadapt)).
	- -> Perform quality filtering ([fastplong](https://github.com/OpenGene/fastplong)).
- **Output**: One cleaned FASTQ file per sample.

### 2. Targeted Read Filtering
- **Input**: The cleaned FASTQ file.
- **Process**:
	- -> Classify reads against a protein database (e.g., all viral proteins) using [DIAMOND](https://github.com/bbuchfink/diamond) BLASTX.
	- -> Extract only the reads that hit the database.
- **Output**: A smaller FASTQ file containing only the reads of interest ("target reads").

### 3. Denovo Assembly
- **Input**: The "target reads" FASTQ file.
- **Process**: The same set of target reads is given to each active assembler.
- **Output**: Multiple raw assembly files (one FASTA per assembler).

### 4. Contig Annotation
- **Input**: The raw assembly files from each assembler.
- **Process**:
	- -> Contigs from all assemblers are collected and annotated against the protein database (DIAMOND).
	- -> A custom script post-processes the results to assign the best taxonomic lineage to each contig.
- **Output**: An annotated table (TSV) for each assembly, linking contigs to potential viruses.

### 5. Per-Assembler Evaluation
- (These steps run independently for each assembler's output to assess its individual quality)
- **Process**:
	- -> General Stats: Calculate basic assembly metrics (N50, total length, etc.).
	- -> Read Mapping: Map all high-quality reads back to the assembly to check coverage.
	- -> Viral QC: Run CheckV to estimate viral genome completeness and check for contamination.
- Output: Individual quality reports and statistics for each assembly.

### 6. Targeted Comparative Analysis
- (This final stage brings the results together for a direct comparison on specific viruses of interest)
- **Process**:
	- -> Binning: Based on the annotations, contigs are binned into separate FASTA files by virus species.
	- -> Comparative QUAST: For each virus of interest (e.g., Monkeypox virus), the binned contigs from all assemblers are compared against the official reference genome.
	- -> Contig Mapping: The binned contigs are also mapped to the reference genome to visualize alignment.
- **Output**: A final, comparative QUAST report showing which assembler performed best for each targeted virus.

# Databases

LoViMAB requires a DIAMOND and CheckV database. Instructions on how to acquire/build these databases:

- [https://bitbucket.org/berkeleylab/checkv](https://bitbucket.org/berkeleylab/checkv)
- [https://github.com/bbuchfink/diamond/wiki](https://github.com/bbuchfink/diamond/wiki)
