--- Version History ---

v1.0.0:
  Initial open-source release of LoViMAB. Workflow includes:
  1. Merge and decompress
  2. Quality control
  3. Pre-filter reads /w diamond blastx database
  4. Run assemblies
  5. Annotate contigs /w diamond blastx database
  6. Post process annotations
  7. Calculate contig stats
  8. Map all QC reads back to assemblies
  9. Run CheckV
  10. Bin contigs by virus (species name) of interest
  11. Map contigs back to viruses of interest
  12. Targeted quast for each virus of interest
