# File: scripts/summarize_checkv.py

import sys
import pandas as pd
from pathlib import Path

def main(checkv_output_path, miuvig_output_path, file_paths):
    """
    Parses multiple CheckV quality_summary.tsv files, creates summary tables,
    and saves them to the specified output paths.
    """

    all_checkv_data = []
    all_miuvig_data = []

    for f_path in file_paths:
        filepath = Path(f_path)
        try:
            # Assumes parent dir is like 'Sample1_metaflye'
            assembler = filepath.parent.name.split('_')[-1]
            df = pd.read_csv(filepath, sep='\t')

            if not df.empty:
                checkv_counts = df['checkv_quality'].value_counts().to_dict()
                checkv_counts['Assembler'] = assembler
                all_checkv_data.append(checkv_counts)

                miuvig_counts = df['miuvig_quality'].value_counts().to_dict()
                miuvig_counts['Assembler'] = assembler
                all_miuvig_data.append(miuvig_counts)
        except Exception as e:
            print(f"Warning: Could not process CheckV file {filepath}. Error: {e}", file=sys.stderr)

    # --- Process and save CheckV quality summary ---
    if all_checkv_data:
        checkv_summary_df = pd.DataFrame(all_checkv_data)
        # Group by assembler and sum the counts across all samples
        checkv_summary_df = checkv_summary_df.groupby('Assembler').sum().reset_index()
        cols = ['Assembler'] + sorted([col for col in checkv_summary_df.columns if col != 'Assembler'])
        checkv_summary_df = checkv_summary_df[cols].fillna(0).astype({col: int for col in cols if col != 'Assembler'})
    else:
        checkv_summary_df = pd.DataFrame(columns=['Assembler'])
    
    # Save to the specific path provided by Snakemake
    checkv_summary_df.to_csv(checkv_output_path, index=False)

    # --- Process and save MIUViG quality summary ---
    if all_miuvig_data:
        miuvig_summary_df = pd.DataFrame(all_miuvig_data)
        # Group by assembler and sum the counts across all samples
        miuvig_summary_df = miuvig_summary_df.groupby('Assembler').sum().reset_index()
        cols = ['Assembler'] + sorted([col for col in miuvig_summary_df.columns if col != 'Assembler'])
        miuvig_summary_df = miuvig_summary_df[cols].fillna(0).astype({col: int for col in cols if col != 'Assembler'})
    else:
        miuvig_summary_df = pd.DataFrame(columns=['Assembler'])
        
    # Save to the specific path provided by Snakemake
    miuvig_summary_df.to_csv(miuvig_output_path, index=False)

if __name__ == "__main__":
    # The script now expects:
    # 1. Path for checkv_summary.csv
    # 2. Path for miuvig_summary.csv
    # 3. All the input .tsv files that follow
    checkv_out = sys.argv[1]
    miuvig_out = sys.argv[2]
    input_files = sys.argv[3:]
    main(checkv_out, miuvig_out, input_files)
