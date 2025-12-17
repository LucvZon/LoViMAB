rule generate_quarto_report:
    input:
        # --- Global Summary Files (for index.qmd) ---
        global_benchmark_summary=expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "benchmark_summary.csv"),assembly_type=ASSEMBLY_TYPES),
        global_assembly_summary=expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "assembly_summary.csv"),assembly_type=ASSEMBLY_TYPES),
        global_checkv_summary=expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "checkv_global_summary.csv"),assembly_type=ASSEMBLY_TYPES),
        global_miuvig_summary=expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "miuvig_global_summary.csv"),assembly_type=ASSEMBLY_TYPES),

        # --- Per-Sample Stat Files (for sample_*.qmd chapters) ---
        per_sample_benchmarks=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_benchmark_summary.csv"), sample=SAMPLES, assembly_type=ASSEMBLY_TYPES),
        per_sample_assemblies=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_assembly_summary.csv"), sample=SAMPLES, assembly_type=ASSEMBLY_TYPES),
        per_sample_checkv_summaries=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_checkv_summary.csv"), sample=SAMPLES, assembly_type=ASSEMBLY_TYPES),
        per_sample_miuvig_summaries=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_miuvig_summary.csv"), sample=SAMPLES, assembly_type=ASSEMBLY_TYPES),
        per_sample_read_pcts=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_reads_leftover_pct.txt"), sample=SAMPLES),
        per_sample_mean_lengths=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_mean_read_length.txt"), sample=SAMPLES),

        # --- Raw Data Files (for on-the-fly calculations in the run block) ---
        qc_reads=expand(os.path.join(QC_DIR, "{sample}.qc.fastq"), sample=SAMPLES),
        target_reads=expand(os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq"), sample=SAMPLES),

        # --- Checkpoint / Flag Files (to ensure upstream rules are complete) ---
        targeted_comparisons_done=os.path.join(STATS_DIR, "all_targeted_comparisons.done"),

        # --- Configuration & Metadata Files ---
        config_file="config/lovimab.yaml",
        software_versions=os.path.join(RESULTS_DIR, "report", "versions.tsv")
    output:
        os.path.join(RESULTS_DIR, "report", "final_summary_report", "index.html")
    params:
        temp_quarto_src=os.path.join(RESULTS_DIR, "report", "quarto_book_src")
    log:
        os.path.join(LOG_DIR, "generate_quarto_report.log")
    run:
        import glob
        import shutil
        from pathlib import Path
        import re
        import gzip
        import pandas as pd # This can be removed I think... 

        # --- Helper function to count reads ---
        def count_fastq_reads(filepath):
            """
            Efficiently counts reads in a FASTQ file, handling both .gz and plain text.
            Reads the file in chunks to be memory-efficient and fast.
            """
            # Choose the correct open function based on file extension
            _open = gzip.open if str(filepath).endswith(".gz") else open
            
            # Open in binary mode ('rb') for speed
            with _open(filepath, "rb") as f:
                # Create a generator to read the file in 1MB chunks
                buf_gen = iter(lambda: f.read(1024 * 1024), b"")
                # Sum the count of newline characters in each chunk
                line_count = sum(buf.count(b"\n") for buf in buf_gen)
            
            # Each FASTQ record has 4 lines
            return line_count // 4

        # --- Helper to sanitize names for IDs/filenames ---
        def sanitize(name):
            return re.sub(r'[^a-zA-Z0-9_]', '_', name)

        # --- Helper: Merge CSVs from different types into one ---
        def merge_and_copy_csvs(file_list, output_filename, dest_dir):
            """
            Reads a list of CSVs, adds an 'Assembly Type' column based on the 
            parent directory name, merges them, and saves to dest_dir.
            """
            dfs = []
            for f in file_list:
                try:
                    # Infer type from parent directory (primary/secondary/final)
                    # structure is .../report/{type}/file.csv or .../per_sample/{type}/file.csv
                    a_type = os.path.basename(os.path.dirname(f))
                    
                    # Read CSV
                    df = pd.read_csv(f)
                    if not df.empty:
                        # Insert Type column at the beginning
                        df.insert(0, "Assembly Type", a_type)
                        dfs.append(df)
                except Exception as e:
                    print(f"Warning: Could not process {f}: {e}")

            dest_path = os.path.join(dest_dir, output_filename)
            if dfs:
                combined_df = pd.concat(dfs, ignore_index=True)
                combined_df.to_csv(dest_path, index=False)
            else:
                # Create empty file with headers if possible, or just touch it
                # Taking a guess at headers based on filename is hard, so we assume
                # the R code handles empty files or we just create an empty CSV.
                with open(dest_path, 'w') as f:
                    f.write("Assembly Type,Process,Note\n") # Minimal header

        # --- Pre-load templates into memory ---
        with open("templates/sample_chapter_template.qmd", 'r') as f:
            main_template_str = f.read()
        # Load the new child template containing the R code
        with open("templates/quast_child_template.qmd", 'r') as f:
            child_template_str = f.read()
        
        # --- Setup Temp Directory ---
        if os.path.exists(params.temp_quarto_src):
            shutil.rmtree(params.temp_quarto_src)
        os.makedirs(params.temp_quarto_src, exist_ok=True)

        # --- Copy static assets ---
        shutil.copy("templates/index.qmd", params.temp_quarto_src)
        shutil.copy("templates/config_chapter.qmd", params.temp_quarto_src)
        shutil.copy(input.software_versions, params.temp_quarto_src)
        shutil.copy(input.config_file, params.temp_quarto_src)

        # Merge Global Summaries
        # This creates a single "benchmark_summary.csv" in the temp dir containing rows for primary, secondary, and final
        merge_and_copy_csvs(input.global_benchmark_summary, "benchmark_summary.csv", params.temp_quarto_src)
        merge_and_copy_csvs(input.global_assembly_summary, "assembly_summary.csv", params.temp_quarto_src)
        merge_and_copy_csvs(input.global_checkv_summary, "checkv_global_summary.csv", params.temp_quarto_src)
        merge_and_copy_csvs(input.global_miuvig_summary, "miuvig_global_summary.csv", params.temp_quarto_src)

        for f in input.per_sample_checkv_summaries: shutil.copy(f, params.temp_quarto_src)
        for f in input.per_sample_miuvig_summaries: shutil.copy(f, params.temp_quarto_src)
        for f in input.per_sample_benchmarks: shutil.copy(f, params.temp_quarto_src)
        for f in input.per_sample_assemblies: shutil.copy(f, params.temp_quarto_src)

        chapter_files = ["index.qmd"]

        # --- Generate Sample Chapters ---
        for i, sample in enumerate(SAMPLES):
            chapter_filename = f"sample_{sample}.qmd"
            chapter_files.append(chapter_filename)
            
            # 1. Gather basic stats
            with open(input.per_sample_read_pcts[i], 'r') as f: reads_leftover_pct = f.read().strip()
            with open(input.per_sample_mean_lengths[i], 'r') as f: mean_read_length = f.read().strip()

            total_read_count = count_fastq_reads(input.qc_reads[i])
            target_read_count = count_fastq_reads(input.target_reads[i])
            
            # 2. Format main template
            chapter_content = main_template_str.format(
                sample_name=sample,
                total_read_count=f"{total_read_count:,}",
                target_read_count=f"{target_read_count:,}",
                reads_leftover_pct=reads_leftover_pct,
                mean_read_length=mean_read_length,
                checkv_csv_filename=os.path.basename(input.per_sample_checkv_summaries[i]),
                miuvig_csv_filename=os.path.basename(input.per_sample_miuvig_summaries[i])
            )

            # 3. Find and Process Dynamic QUAST Results
            # Check if we have ANY reports to decide if we print the main header
            any_reports_exist = False
            for a_type in ASSEMBLY_TYPES:
                 if glob.glob(os.path.join(STATS_DIR, "targeted_quast", a_type, sample, "*", "report.tsv")):
                     any_reports_exist = True
                     break
            
            if any_reports_exist:
                chapter_content += "\n\n## Targeted Assembly Quality (QUAST)\n"

                # Loop through types to keep order: Primary -> Secondary -> Final
                for a_type in ASSEMBLY_TYPES:
                    
                    # Pattern matches reports specifically for THIS assembly type
                    type_pattern = os.path.join(STATS_DIR, "targeted_quast", a_type, sample, "*", "report.tsv")
                    found_reports = sorted(glob.glob(type_pattern))

                    if not found_reports:
                        continue
                        
                    # Create the title string (e.g., "Primary Assembly")
                    assembly_title = f"{a_type.capitalize()} Assembly"

                    for report_path_str in found_reports:
                        report_path = Path(report_path_str)
                        # Structure: .../{a_type}/{sample}/{virus_folder}/report.tsv
                        virus_folder_name = report_path.parent.name
                        virus_display = virus_folder_name.replace("_", " ")

                        # Create unique filename including assembly type to avoid collisions
                        # e.g., Sample1_primary_Human_alpha_report.tsv
                        unique_report_filename = f"{sample}_{a_type}_{virus_folder_name}_report.tsv"
                        dest_path = os.path.join(params.temp_quarto_src, unique_report_filename)
                        
                        shutil.copy(report_path, dest_path)

                        # Unique ID for R chunk labels
                        unique_id = sanitize(f"{sample}_{a_type}_{virus_folder_name}")

                        # Format template
                        formatted_child = child_template_str.format(
                            virus_display_name=virus_display,
                            unique_id=unique_id,
                            report_filename=unique_report_filename,
                            assembly_type_name=assembly_title
                        )
                        chapter_content += "\n" + formatted_child + "\n"

            # 4. Write the final generated chapter
            with open(os.path.join(params.temp_quarto_src, chapter_filename), 'w') as f:
                f.write(chapter_content)
				
        chapter_files.append("config_chapter.qmd")

        # --- Generate _quarto.yml ---
        yaml_content = """project:
  type: book
  output-dir: _book
book:
  title: "LoViMAB Workflow Summary"
  favicon: images/favicon.png
  chapters:
"""
        for chapter in chapter_files: yaml_content += f"    - {chapter}\n"
        yaml_content += """
format:
  html:
    grid:
      sidebar-width: 200px
      body-width: 1250px
      margin-width: 200px
      gutter-width: 1em
    theme: cyborg
    toc: true
engine: knitr
"""
        with open(os.path.join(params.temp_quarto_src, "_quarto.yml"), 'w') as f: f.write(yaml_content)
        
        # Add simple CSS for table scrolling if they get wide
        with open(os.path.join(params.temp_quarto_src, "styles.css"), 'w') as f:
            f.write(".cell-output-display { overflow-x: auto; }")

        # --- Render ---
        # Use absolute path for log to avoid 'cd' issues
        abs_log_path = os.path.abspath(str(log))
        shell(f"quarto render {params.temp_quarto_src} --to html &> {abs_log_path}")
        
        # --- Finalize ---
        final_output_dir = os.path.dirname(output[0])
        os.makedirs(final_output_dir, exist_ok=True)
        
        # Use python to move files
        book_output = os.path.join(params.temp_quarto_src, "_book")
        for item in os.listdir(book_output):
            s = os.path.join(book_output, item)
            d = os.path.join(final_output_dir, item)
            if os.path.exists(d):
                if os.path.isdir(d): shutil.rmtree(d)
                else: os.remove(d)
            shutil.move(s, d)
            
        shutil.rmtree(params.temp_quarto_src)
