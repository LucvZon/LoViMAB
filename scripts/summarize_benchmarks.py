# File: scripts/summarize_benchmarks.py

import sys
import pandas as pd
from pathlib import Path
import subprocess

def process_benchmark_log(filepath_str):
    """Parses a single Snakemake benchmark log file."""
    filepath = Path(filepath_str)
    COLUMNS_TO_PARSE = {'s': float, 'max_rss': float, 'mean_load': float}
    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=0,
            usecols=COLUMNS_TO_PARSE.keys(),
            dtype=COLUMNS_TO_PARSE,
            comment='#'
        )
        if not df.empty:
            stats = df.iloc[0].to_dict()
            stats['process'] = filepath.parent.name
            return stats
    except Exception as e:
        print(f"Warning: Could not parse benchmark file {filepath}. Error: {e}", file=sys.stderr)
    return None

def process_stats_file(filepath_str):
    """Parses a single stats.sh output file."""
    filepath = Path(filepath_str)
    try:
        # Read the tab-separated file
        df = pd.read_csv(filepath, sep='\t', header=0)
        
        # 1. Get the stem: "Sample1_metaflye.stats"
        # 2. Remove the unwanted ".stats" suffix
        # 3. Then, get the assembler name
        clean_stem = filepath.stem.removesuffix('.stats') # This is a Python 3.9+ feature
        assembler_name = clean_stem.split('_')[-1]
        
        # Get the first (and only) row of data
        stats = df.iloc[0].to_dict()
        
        # Return a dictionary with the process name and relevant stats
        return {
            'process': assembler_name,
            '# Contigs': stats['n_contigs'],
            'Total Length (bp)': stats['contig_bp'],
            'Largest Contig (bp)': stats['ctg_max']
        }
    except Exception as e:
        print(f"Warning: Could not parse stats file {filepath}. Error: {e}", file=sys.stderr)
    return None

def process_bam_file(filepath_str):
    """Gets total and mapped read counts from a BAM file using samtools."""
    filepath = Path(filepath_str)
    try:
        # Filename format is now "{sample}_{assembler}.bam", so this logic is slightly simpler
        assembler = filepath.stem.split('_')[-1]
        total_cmd = f"samtools view -c {filepath}"
        mapped_cmd = f"samtools view -c -F 4 {filepath}" # -F 4 means "not unmapped"
        
        total_reads = int(subprocess.check_output(total_cmd, shell=True).strip())
        mapped_reads = int(subprocess.check_output(mapped_cmd, shell=True).strip())
        
        return {'process': assembler, 'total_reads': total_reads, 'mapped_reads': mapped_reads}
    except Exception as e:
        print(f"Warning: Could not process BAM file {filepath}. Error: {e}", file=sys.stderr)
    return None

def main(all_files):
    benchmark_data, stats_data, bam_data = [], [], []

    for f in all_files:
        if 'benchmarks' in f and f.endswith('.log'):
            res = process_benchmark_log(f)
            if res: benchmark_data.append(res)
        elif f.endswith('.stats.txt'):
            res = process_stats_file(f)
            if res: stats_data.append(res)
        elif f.endswith('.bam'):
            res = process_bam_file(f)
            if res: bam_data.append(res)
            
    # --- Process Benchmarks ---
    if benchmark_data:
        benchmark_df = pd.DataFrame(benchmark_data)
        summary_df = benchmark_df.groupby('process').mean().reset_index()
        summary_df['Run Time (minutes)'] = (summary_df['s'] / 60).round(2)
        summary_df['Peak Memory (GB)'] = (summary_df['max_rss'] / 1024).round(2)
        summary_df['Mean CPU Cores'] = (summary_df['mean_load'] / 100).round(2)
        benchmark_summary = summary_df[['process', 'Run Time (minutes)', 'Peak Memory (GB)', 'Mean CPU Cores']]
    else:
        benchmark_summary = pd.DataFrame(columns=['process', 'Run Time (minutes)', 'Peak Memory (GB)', 'Mean CPU Cores'])
        
    benchmark_summary.to_csv("benchmark_summary.csv", index=False)

    # --- Process Assembly Outcomes ---
    if stats_data:
        stats_df = pd.DataFrame(stats_data)
        # Average the stats across all samples for each assembler
        assembly_summary = stats_df.groupby('process').mean().reset_index()
        
        # Merge the read mapping statistics.
        if bam_data:
            bam_df = pd.DataFrame(bam_data)
            # For reads, we want the SUM across all samples, not the mean.
            bam_summary = bam_df.groupby('process').sum(numeric_only=True).reset_index()
            bam_summary['Reads Mapped (%)'] = ((bam_summary['mapped_reads'] / bam_summary['total_reads']) * 100).round(2)
            
            assembly_summary = pd.merge(assembly_summary, bam_summary[['process', 'Reads Mapped (%)']], on='process', how='left')
        else:
            assembly_summary['Reads Mapped (%)'] = 0.0 # Or None
    else:
        # Create an empty dataframe if no stats files were found
        assembly_summary = pd.DataFrame(columns=['process', '# Contigs', 'Largest Contig (bp)', 'Total Length (bp)', 'Reads Mapped (%)'])
    
    # Ensure all expected columns are present for consistent output
    final_cols = ['process', '# Contigs', 'Largest Contig (bp)', 'Total Length (bp)', 'Reads Mapped (%)']
    for col in final_cols:
        if col not in assembly_summary.columns:
            assembly_summary[col] = 0.0 # or None
    
    assembly_summary.to_csv("assembly_summary.csv", index=False)

if __name__ == "__main__":
    input_files = sys.argv[1:]
    main(input_files)
