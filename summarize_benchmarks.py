# File: scripts/summarize_benchmarks.py

import sys
import pandas as pd
from pathlib import Path
import subprocess

def process_benchmark_log(filepath_str):
    """Parses a single Snakemake benchmark log file."""
    filepath = Path(filepath_str) # Convert string to Path object
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
            # The process name is derived from the *relative* path that Snakemake knows.
            # We assume the path looks like 'benchmarks/rule_name/sample.log'
            stats['process'] = filepath.parent.name
            return stats
    except Exception as e:
        print(f"Warning: Could not parse benchmark file {filepath}. Error: {e}", file=sys.stderr)
    return None

def process_quast_report(filepath_str):
    """Parses a single QUAST transposed_report.tsv file."""
    filepath = Path(filepath_str)
    try:
        df = pd.read_csv(filepath, sep='\t', header=0)
        df = df.rename(columns={'Assembly': 'process'})
        # Ensure only relevant columns are selected and return as list of dicts
        return df[['process', '# contigs', 'Largest contig', 'Total length']].to_dict('records')
    except Exception as e:
        print(f"Warning: Could not parse QUAST report {filepath}. Error: {e}", file=sys.stderr)
    return []

def process_bam_file(filepath_str):
    """Gets total and mapped read counts from a BAM file using samtools."""
    filepath = Path(filepath_str)
    try:
        assembler = filepath.stem.split('_')[-1]
        total_cmd = f"samtools view -c {filepath}"
        mapped_cmd = f"samtools view -c -F 4 {filepath}" # -F 4 means "not unmapped"
        
        # Use shell=True for simpler command execution, but ensure subprocess security if input is untrusted.
        # Here, paths are from trusted Snakemake, so it's fine.
        total_reads = int(subprocess.check_output(total_cmd, shell=True).strip())
        mapped_reads = int(subprocess.check_output(mapped_cmd, shell=True).strip())
        return {'process': assembler, 'total_reads': total_reads, 'mapped_reads': mapped_reads}
    except Exception as e:
        print(f"Warning: Could not process BAM file {filepath}. Error: {e}", file=sys.stderr)
    return None

def main(all_files):
    benchmark_data, quast_data_list, bam_data = [], [], []

    for f in all_files:
        if 'benchmarks' in f and f.endswith('.log'):
            res = process_benchmark_log(f)
            if res: benchmark_data.append(res)
        elif 'transposed_report.tsv' in f:
            quast_data_list.extend(process_quast_report(f))
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
    if quast_data_list:
        quast_df = pd.DataFrame(quast_data_list)
        assembly_summary = quast_df.groupby('process').mean().reset_index() # Average across samples/viruses
        
        if bam_data:
            bam_df = pd.DataFrame(bam_data)
            bam_summary = bam_df.groupby('process').sum(numeric_only=True).reset_index()
            bam_summary['Reads Mapped (%)'] = (bam_summary['mapped_reads'] / bam_summary['total_reads'] * 100).round(2)
            assembly_summary = pd.merge(assembly_summary, bam_summary[['process', 'Reads Mapped (%)']], on='process', how='left')
        else:
            assembly_summary['Reads Mapped (%)'] = None # Or 0, depending on desired empty state
    else:
        assembly_summary = pd.DataFrame(columns=['process', '# Contigs', 'Largest Contig (bp)', 'Total Length (bp)', 'Reads Mapped (%)'])

    assembly_summary = assembly_summary.rename(columns={'# contigs': '# Contigs', 'Largest contig': 'Largest Contig (bp)', 'Total length': 'Total Length (bp)'})
    
    # Ensure all expected columns are present, even if no data, for consistent CSV output
    final_cols = ['process', '# Contigs', 'Largest Contig (bp)', 'Total Length (bp)', 'Reads Mapped (%)']
    for col in final_cols:
        if col not in assembly_summary.columns:
            assembly_summary[col] = None # or 0, depending on what's desired for empty cells
    
    assembly_summary.to_csv("assembly_summary.csv", index=False)

if __name__ == "__main__":
    input_files = sys.argv[1:]
    main(input_files)
