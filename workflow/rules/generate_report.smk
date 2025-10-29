rule generate_quarto_report:
    input:
        # --- Global Summary Files (for index.qmd) ---
        global_benchmark_summary=os.path.join(RESULTS_DIR, "report", "benchmark_summary.csv"),
        global_assembly_summary=os.path.join(RESULTS_DIR, "report", "assembly_summary.csv"),
        global_checkv_summary=os.path.join(RESULTS_DIR, "report", "checkv_global_summary.csv"),
        global_miuvig_summary=os.path.join(RESULTS_DIR, "report", "miuvig_global_summary.csv"),

        # --- Per-Sample Stat Files (for sample_*.qmd chapters) ---
        per_sample_benchmarks=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_benchmark_summary.csv"), sample=SAMPLES),
        per_sample_assemblies=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_assembly_summary.csv"), sample=SAMPLES),
        per_sample_checkv_summaries=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_checkv_summary.csv"), sample=SAMPLES),
        per_sample_miuvig_summaries=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_miuvig_summary.csv"), sample=SAMPLES),
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
        import gzip # Import the gzip library

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
        shutil.copy(input.global_benchmark_summary, params.temp_quarto_src)
        shutil.copy(input.global_assembly_summary, params.temp_quarto_src)
        shutil.copy(input.global_checkv_summary, params.temp_quarto_src)
        shutil.copy(input.global_miuvig_summary, params.temp_quarto_src)
        shutil.copy(input.software_versions, params.temp_quarto_src)
        shutil.copy(input.config_file, params.temp_quarto_src)
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
            quast_report_pattern = os.path.join(STATS_DIR, "targeted_quast", sample, "*", "report.tsv")
            found_reports = sorted(glob.glob(quast_report_pattern))

            if found_reports:
                chapter_content += "\n\n## Targeted Assembly Quality (QUAST)\n"

                for report_path_str in found_reports:
                    report_path = Path(report_path_str)
                    # Directory structure is .../{sample}/{virus_folder}/report.tsv
                    virus_folder_name = report_path.parent.name
                    virus_display = virus_folder_name.replace("_", " ")
                    
                    # Create a unique filename for this report in the temp dir
                    # e.g., Sample1_Human_alphaherpesvirus_1_report.tsv
                    unique_report_filename = f"{sample}_{virus_folder_name}_report.tsv"
                    dest_path = os.path.join(params.temp_quarto_src, unique_report_filename)
                    
                    # Copy the data file
                    shutil.copy(report_path, dest_path)

                    # Creating a unique ID for Quarto labels (no spaces allowed)
                    unique_id = sanitize(f"{sample}_{virus_folder_name}")

                    # Format the CHILD template and append it
                    formatted_child = child_template_str.format(
                        virus_display_name=virus_display,
                        unique_id=unique_id,
                        report_filename=unique_report_filename
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
