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

snakemake --snakefile workflow/snakefile.smk
```

