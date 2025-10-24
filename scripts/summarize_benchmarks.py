# File: scripts/summarize_benchmarks.py

import sys
import pandas as pd
from pathlib import Path
import subprocess
import argparse

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
            stats['sample'] = filepath.stem # e.g., "Sample1" from "Sample1.log"
            return stats
    except Exception as e:
        print(f"Warning: Could not parse benchmark file {filepath}. Error: {e}", file=sys.stderr)
    return None

def process_stats_file(filepath_str):
    """Parses a single stats.sh output file."""
    filepath = Path(filepath_str)
    try:
        df = pd.read_csv(filepath, sep='\t', header=0)
        
        # Stem is "Sample1_metaflye.stats", so we remove suffix to get "Sample1_metaflye"
        clean_stem = filepath.stem.removesuffix('.stats')
        assembler_name = clean_stem.split('_')[-1]
        sample_name = clean_stem.removesuffix(f"_{assembler_name}")

        stats = df.iloc[0].to_dict()
        
        return {
            'process': assembler_name,
            'sample': sample_name,
            '# Contigs': stats['n_contigs'],
            'Total Length (bp)': stats['contig_bp'],
            'Largest Contig (bp)': stats['ctg_max']
        }
    except Exception as e:
        print(f"Warning: Could not parse stats file {filepath}. Error: {e}", file=sys.stderr)
    return None

def process_bam_file(filepath_str, threads=1):
    """Gets total and mapped read counts from a BAM file using samtools."""
    filepath = Path(filepath_str)
    try:
        assembler = filepath.stem.split('_')[-1]
        sample = filepath.stem.removesuffix(f"_{assembler}")
        total_cmd = f"samtools view -@ {threads} -c {filepath}"
        mapped_cmd = f"samtools view -@ {threads} -c -F 4 {filepath}"
        
        total_reads = int(subprocess.check_output(total_cmd, shell=True).strip())
        mapped_reads = int(subprocess.check_output(mapped_cmd, shell=True).strip())
        
        return {'process': assembler, 'sample': sample, 'total_reads': total_reads, 'mapped_reads': mapped_reads}
    except Exception as e:
        print(f"Warning: Could not process BAM file {filepath}. Error: {e}", file=sys.stderr)
    return None

def main(all_files, threads):
    benchmark_data, stats_data, bam_data = [], [], []

    for f in all_files:
        if 'benchmarks' in f and f.endswith('.log'):
            res = process_benchmark_log(f)
            if res: benchmark_data.append(res)
        elif f.endswith('.stats.txt'):
            res = process_stats_file(f)
            if res: stats_data.append(res)
        elif f.endswith('.bam'):
            res = process_bam_file(f, threads)
            if res: bam_data.append(res)
            
    # --- Convert to DataFrames ---
    benchmark_df = pd.DataFrame(benchmark_data) if benchmark_data else pd.DataFrame(columns=['sample', 'process', 's', 'max_rss', 'mean_load'])
    stats_df = pd.DataFrame(stats_data) if stats_data else pd.DataFrame(columns=['sample', 'process', '# Contigs', 'Total Length (bp)', 'Largest Contig (bp)'])
    bam_df = pd.DataFrame(bam_data) if bam_data else pd.DataFrame(columns=['sample', 'process', 'total_reads', 'mapped_reads'])

    # --- Process Overall (Averaged) Summaries ---
    if not benchmark_df.empty:
        avg_benchmark_df = benchmark_df.groupby('process').mean(numeric_only=True).reset_index()
        avg_benchmark_df['Run Time (minutes)'] = (avg_benchmark_df['s'] / 60).round(2)
        avg_benchmark_df['Peak Memory (GB)'] = (avg_benchmark_df['max_rss'] / 1024).round(2)
        avg_benchmark_df['Mean CPU Cores'] = (avg_benchmark_df['mean_load'] / 100).round(2)
        benchmark_summary = avg_benchmark_df[['process', 'Run Time (minutes)', 'Peak Memory (GB)', 'Mean CPU Cores']]
    else:
        benchmark_summary = pd.DataFrame(columns=['process', 'Run Time (minutes)', 'Peak Memory (GB)', 'Mean CPU Cores'])
    benchmark_summary.to_csv("benchmark_summary.csv", index=False)

    assembly_summary_avg = pd.DataFrame()
    if not stats_df.empty:
        assembly_summary_avg = stats_df.groupby('process').mean(numeric_only=True).reset_index()
        if not bam_df.empty:
            sum_bam_df = bam_df.groupby('process').sum(numeric_only=True).reset_index()
            if 'total_reads' in sum_bam_df.columns and sum_bam_df['total_reads'].sum() > 0:
                sum_bam_df['Reads Mapped (%)'] = ((sum_bam_df['mapped_reads'] / sum_bam_df['total_reads']) * 100).round(2)
                assembly_summary_avg = pd.merge(assembly_summary_avg, sum_bam_df[['process', 'Reads Mapped (%)']], on='process', how='left')
    
    final_cols = ['process', '# Contigs', 'Largest Contig (bp)', 'Total Length (bp)', 'Reads Mapped (%)']
    for col in final_cols:
        if col not in assembly_summary_avg.columns:
            assembly_summary_avg[col] = 0.0
    assembly_summary_avg[final_cols].to_csv("assembly_summary.csv", index=False)

    # --- Process Per-Sample Summaries ---
    all_samples = pd.concat([benchmark_df.get('sample'), stats_df.get('sample'), bam_df.get('sample')]).dropna().unique()

    for sample in all_samples:
        # Filter data for the current sample
        sample_benchmark_df = benchmark_df[benchmark_df['sample'] == sample].copy()
        sample_stats_df = stats_df[stats_df['sample'] == sample].copy()
        sample_bam_df = bam_df[bam_df['sample'] == sample].copy()

        # Generate benchmark summary for the sample
        if not sample_benchmark_df.empty:
            sample_benchmark_df['Run Time (minutes)'] = (sample_benchmark_df['s'] / 60).round(2)
            sample_benchmark_df['Peak Memory (GB)'] = (sample_benchmark_df['max_rss'] / 1024).round(2)
            sample_benchmark_df['Mean CPU Cores'] = (sample_benchmark_df['mean_load'] / 100).round(2)
            sample_benchmark_summary = sample_benchmark_df[['process', 'Run Time (minutes)', 'Peak Memory (GB)', 'Mean CPU Cores']]
            sample_benchmark_summary.to_csv(f"{sample}_benchmark_summary.csv", index=False)

        # Generate assembly summary for the sample
        if not sample_stats_df.empty:
            assembly_summary_persample = sample_stats_df.copy()
            if not sample_bam_df.empty and 'total_reads' in sample_bam_df.columns and sample_bam_df['total_reads'].sum() > 0:
                sample_bam_df['Reads Mapped (%)'] = ((sample_bam_df['mapped_reads'] / sample_bam_df['total_reads']) * 100).round(2)
                assembly_summary_persample = pd.merge(assembly_summary_persample, sample_bam_df[['process', 'Reads Mapped (%)']], on='process', how='left')
            else:
                assembly_summary_persample['Reads Mapped (%)'] = 0.0
            
            for col in final_cols:
                if col not in assembly_summary_persample.columns:
                    assembly_summary_persample[col] = 0.0
            assembly_summary_persample[final_cols].to_csv(f"{sample}_assembly_summary.csv", index=False)


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Summarize benchmark and assembly results.")
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for samtools.')
    parser.add_argument('input_files', nargs='+', help='List of input files (benchmarks, stats, bams).')
    
    args = parser.parse_args()
    
    # Pass arguments to the main function
    main(args.input_files, args.threads)
