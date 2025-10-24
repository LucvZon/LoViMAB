# ===================================================================
#          SNAKEMAKE WORKFLOW FOR TARGETED METAGENOMIC ASSEMBLY
# ===================================================================

# VERSION: 1.0.0

import os
import glob

# --- Load Configuration ---
configfile: "config/lovimab.yaml"

# --- Onstart: Validate Configuration ---
onstart:
    # Check if every file path listed in configfile exists.
    for virus_name, ref_path in config["viruses_of_interest"].items():
        if not os.path.exists(ref_path):
            raise WorkflowError(
                f"Configuration Error: The reference genome path for virus '{virus_name}' "
                f"does not exist. Please check the path in your config file.\n"
                f"  -> Path provided: {ref_path}"
            )
    for database, db_path in config["paths"].items():
        if not os.path.exists(db_path):
            raise WorkflowError(
                f"Configuration Error: The following database '{db_path}' "
                f"does not exist. Please check the path in your config file.\n"
                f"  -> Path provided: {db_path}"
            )

# --- Global Variables ---
SAMPLES = list(config["samples"].keys())
ASSEMBLERS_CONFIG = config["assemblers"]
ACTIVE_ASSEMBLERS = [asm for asm, active in ASSEMBLERS_CONFIG.items() if active]

# Define output directories
RESULTS_DIR = "results"
QC_DIR = os.path.join(RESULTS_DIR, "1_quality_control")
READ_CLASSIFICATION_DIR = os.path.join(RESULTS_DIR, "2_read_classification")
ASSEMBLY_DIR = os.path.join(RESULTS_DIR, "3_assemblies")
ANNOTATION_DIR = os.path.join(RESULTS_DIR, "4_annotation")
STATS_DIR = os.path.join(RESULTS_DIR, "5_stats_and_qc")
LOG_DIR = "logs"
BENCH_DIR = "benchmarks"

# --- Wildcard Constraints ---
# Prevent ambiguous matching by telling Snakemake exactly what the
# 'assembler' wildcard can be.
wildcard_constraints:
    assembler="|".join(ACTIVE_ASSEMBLERS)

# --- Helper function to get correct assembly output file ---
def get_assembly_fasta(wildcards):
    assembler = wildcards.assembler
    if assembler == "metaflye":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "penguin":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "raven":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "canu":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "canu_assembly.contigs.fasta")
    elif assembler == "myloasm":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "assembly_primary.fa")
    elif assembler == "metamdbg":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "wtdbg2":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "shasta":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "Assembly.fasta")
    elif assembler == "miniasm":
        return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "final_assembly.fasta")
    # Add other assemblers here if needed in the future

# ===================================================================
#                              TARGET RULE
# ===================================================================
rule all:
    input:
        # --- Generic results ---
        expand(os.path.join(STATS_DIR, "contig_stats", "{sample}_{assembler}.stats.txt"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        expand(os.path.join(STATS_DIR, "reads_to_contigs", "{sample}_{assembler}.bam"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        expand(os.path.join(STATS_DIR, "checkv", "{sample}_{assembler}"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        expand(os.path.join(ANNOTATION_DIR, "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        # --- Specific, Targeted Reports ---
        os.path.join(STATS_DIR, "all_targeted_comparisons.done"),
        os.path.join(STATS_DIR, "all_contig_mappings.done"),
        os.path.join(RESULTS_DIR, "report", "final_summary_report", "index.html")


# ===================================================================
#                          WORKFLOW RULES
# ===================================================================

def get_possible_quast_comparisons(wildcards):
    """
    This function correctly scans the output of the 'bin_contigs_by_taxonomy' checkpoint
    to determine which grouped QUAST comparison reports are possible to generate.
    """
    found_combos = set()
    
    # Iterate through each sample and assembler combination
    # and get the checkpoint output for each one individually.
    for sample in SAMPLES:
        for assembler in ACTIVE_ASSEMBLERS:
            # Get the output directory from the checkpoint for a single S/A combo.
            binned_dir = checkpoints.bin_contigs_by_taxonomy.get(
                sample=sample, assembler=assembler
            ).output.binned_dir
            
            # Scan this specific directory for any .fasta files it created.
            for fasta_file in glob.glob(os.path.join(str(binned_dir), "*.fasta")):
                virus_name = os.path.basename(fasta_file).replace(".fasta", "")
                
                # Check if this virus is one we are interested in.
                if virus_name in config["viruses_of_interest"]:
                    # Store the unique (sample, virus) combination.
                    found_combos.add((sample, virus_name))

    # Build the final list of report paths from the unique combos found.
    possible_reports = []
    for sample, virus_name in found_combos:
        possible_reports.append(
            os.path.join(STATS_DIR, "targeted_quast", sample, virus_name, "report.html")
        )
        
    return possible_reports


def get_possible_contig_mappings(wildcards):
    """
    This function scans the checkpoint output to determine which binned contigs
    can be mapped to their corresponding reference of interest.
    """
    possible_bam_files = []
    
    for sample in SAMPLES:
        for assembler in ACTIVE_ASSEMBLERS:
            binned_dir = checkpoints.bin_contigs_by_taxonomy.get(
                sample=sample, assembler=assembler
            ).output.binned_dir
            
            for fasta_file in glob.glob(os.path.join(str(binned_dir), "*.fasta")):
                virus_name = os.path.basename(fasta_file).replace(".fasta", "")
                
                # Check if this discovered virus is in the high-priority list.
                if virus_name in config["viruses_of_interest"]:
                    possible_bam_files.append(
                        os.path.join(STATS_DIR, "contigs_to_ref", sample, assembler, f"{virus_name}.bam")
                    )
                    
    return possible_bam_files


# Step 1: Merge raw reads for a given sample
rule merge_reads:
    output:
        temp(os.path.join(QC_DIR, "{sample}.merged.fastq"))
    params:
        prefix=lambda wildcards: config["samples"][wildcards.sample]
    log:
        os.path.join(LOG_DIR, "merge_reads", "{sample}.log")
    shell:
        "zcat {params.prefix}*.fastq* > {output} 2> {log}"

# Step ?: Trim adapter sequences 
rule trim_adapters:
    input:
        os.path.join(QC_DIR, "{sample}.merged.fastq")
    output:
        fastq=os.path.join(QC_DIR, "{sample}.trimmed.fastq"),
        tmp_out=temp("{sample}.cutadapt.tmp")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "trimming", "{sample}.log")
    shell:
        """
        cutadapt -j {threads} -e 0.2 -n 5 -m 150 --revcomp -a GTTTCCCACTGGAGGATA...TATCCTCCAGTGGGAAAC {input} > {output.tmp_out} 2> {log}
        cutadapt -j {threads} -u 9 -u -9 {output.tmp_out} > {output.fastq} 2>> {log}
        """

# Step 2: Perform quality control on merged reads
rule quality_control:
    input:
        os.path.join(QC_DIR, "{sample}.trimmed.fastq")
    output:
        fastq=os.path.join(QC_DIR, "{sample}.qc.fastq"),
        json=os.path.join(QC_DIR, "{sample}.qc.json"),
        html=os.path.join(QC_DIR, "{sample}.qc.html")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "quality_control", "{sample}.log")
    shell:
        "fastplong -i {input} -o {output.fastq} "
        "--length_required 150 --qualified_quality_phred 10 -j {output.json} -h {output.html} "
        "--unqualified_percent_limit 50 --disable_adapter_trimming --thread {threads} &> {log}"

# Step 3: Classify reads with DIAMOND against a custom database
rule classify_reads_diamond:
    input:
        os.path.join(QC_DIR, "{sample}.qc.fastq")
    output:
        os.path.join(READ_CLASSIFICATION_DIR, "{sample}.diamond_annotation.tsv")
    params:
        db=config["paths"]["diamond_db"],
        sensitivity_flag=f'--{config["params"]["diamond_sensitivity"]}'
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "classify_reads_diamond", "{sample}.log")
    benchmark:
        os.path.join(BENCH_DIR, "classify_reads_diamond", "{sample}.log")
    shell:
        "diamond blastx {params.sensitivity_flag} -d {params.db} -q {input} -o {output} "
        "-f 6 qseqid -k 1 --threads {threads} &> {log}" # Using -k 1 for best hit

# Step 4: Extract target reads based on DIAMOND classification
rule extract_target_reads:
    input:
        reads=os.path.join(QC_DIR, "{sample}.qc.fastq"),
        ids=os.path.join(READ_CLASSIFICATION_DIR, "{sample}.diamond_annotation.tsv")
    output:
        os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
    log:
        os.path.join(LOG_DIR, "extract_target_reads", "{sample}.log")
    shell:
        "awk '{{print $1}}' {input.ids} | sort -u > {output}.ids.txt 2> {log}; "
        "seqtk subseq {input.reads} {output}.ids.txt > {output} 2>> {log}"
		
# Step 4b: Calculate the percentage of reads kept after classification
rule calculate_reads_leftover:
    input:
        total_reads_file=os.path.join(QC_DIR, "{sample}.qc.fastq"),
        target_reads_file=os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
    output:
        os.path.join(STATS_DIR, "per_sample", "{sample}_reads_leftover_pct.txt")
    log:
        os.path.join(LOG_DIR, "calculate_reads_leftover", "{sample}.log")
    shell:
        """
        # Count lines in each file and divide by 4 to get the number of reads
        total_reads=$( (wc -l < {input.total_reads_file} || true) | awk '{{print $1 / 4}}' )
        target_reads=$( (wc -l < {input.target_reads_file} || true) | awk '{{print $1 / 4}}' )

        # Use awk for safe floating-point arithmetic, even if counts are zero
        awk -v target="$target_reads" -v total="$total_reads" \
        'BEGIN {{ if (total > 0) {{ printf "%.2f", (target / total) * 100 }} else {{ print "0.00" }} }}' > {output} 2> {log}
        """

# --- ASSEMBLY RULES ---
include: "rules/assemblers.smk"

# --- POST-ASSEMBLY PROCESSING ---

# Step 6a: Rename contigs to ensure uniqueness before aggregation
rule rename_contigs:
    input:
        fasta=get_assembly_fasta
    output:
        os.path.join(ASSEMBLY_DIR, "{sample}", "{assembler}", "renamed_contigs.fasta")
    params:
        prefix="{assembler}_"
    shell:
        "seqkit replace -p '^' -r '{params.prefix}' {input.fasta} > {output}"

# Step 6b: Aggregate all renamed contigs from a single sample
rule aggregate_contigs:
    input:
        expand(os.path.join(ASSEMBLY_DIR, "{{sample}}", "{assembler}", "renamed_contigs.fasta"), assembler=ACTIVE_ASSEMBLERS)
    output:
        os.path.join(ANNOTATION_DIR, "{sample}", "aggregated_contigs.fasta")
    shell:
        "cat {input} > {output}"

# Step 6c: Annotate the single aggregated contig file
rule annotate_aggregated_contigs:
    input:
        os.path.join(ANNOTATION_DIR, "{sample}", "aggregated_contigs.fasta")
    output:
        os.path.join(ANNOTATION_DIR, "{sample}", "aggregated_annotation.tsv")
    params:
        db=config["paths"]["diamond_db"]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "annotate_aggregated", "{sample}.log")
    shell:
        "diamond blastx -d {params.db} -q {input} -o {output} "
        "-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids "
        "--threads {threads} -b 10 -c 1 &> {log}"

# Step 7: Split the aggregated annotation file back into per-assembler files
rule split_annotations:
    input:
        os.path.join(ANNOTATION_DIR, "{sample}", "aggregated_annotation.tsv")
    output:
        os.path.join(ANNOTATION_DIR, "{sample}", "split", "{assembler}_annotation.tsv")
    params:
        prefix="{assembler}_"
    shell:
        "grep -E '^{params.prefix}' {input} > {output} || touch {output}" # touch avoids error if no hits found

# Step 8: Post-process the per-assembler annotation file
rule post_process_annotation:
    input:
        annotation=os.path.join(ANNOTATION_DIR, "{sample}", "split", "{assembler}_annotation.tsv"),
        contigs=os.path.join(ASSEMBLY_DIR, "{sample}", "{assembler}", "renamed_contigs.fasta")
    output:
        annotated=os.path.join(ANNOTATION_DIR, "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"),
        unannotated=os.path.join(ANNOTATION_DIR, "{sample}", "post_processed", "{assembler}_unannotated_contigs.tsv")
    params:
        script=workflow.source_path("../scripts/post_process_diamond_v1.0.py")
    log:
        os.path.join(LOG_DIR, "post_process", "{sample}_{assembler}.log")
    shell:
        "python {params.script} -i {input.annotation} -c {input.contigs} "
        "-o {output.annotated} -u {output.unannotated} -log {log}"

# EXPERIMENTAL STEP
checkpoint bin_contigs_by_taxonomy:
    input:
        annotation=os.path.join(ANNOTATION_DIR, "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"),
        contigs=os.path.join(ASSEMBLY_DIR, "{sample}", "{assembler}", "renamed_contigs.fasta")
    output:
        binned_dir=directory(os.path.join(ANNOTATION_DIR, "{sample}", "binned_contigs", "{assembler}"))
    params:
        species_col=14, # Column name containing species information
        contig_col=1
    log:
        os.path.join(LOG_DIR, "bin_contigs", "{sample}_{assembler}.log")
    shell:
        """
        # Ensure the output directory exists
        mkdir -p {output.binned_dir}

        # Read the annotation file line by line (skipping the header)
        tail -n +2 {input.annotation} | while IFS=$'\\t' read -r -a line; do
            contig_id="${{line[{params.contig_col}-1]}}"
            
            # Sanitize the species name to create a valid filename (e.g., "Human alphaherpesvirus 1" -> "Human_alphaherpesvirus_1.fasta")
            species_name=$(echo "${{line[{params.species_col}-1]}}" | sed 's/ /_/g; s/[^a-zA-Z0-9_.-]//g')
            
            # If the species name is not empty, extract the contig and append it to the corresponding species FASTA file
            if [[ -n "$species_name" ]]; then
                seqtk subseq {input.contigs} <(echo "$contig_id") >> {output.binned_dir}/$species_name.fasta
            fi
        done > {log} 2>&1
        """


rule aggregate_targeted_comparisons:
    input:
        get_possible_quast_comparisons
    output:
        touch(os.path.join(STATS_DIR, "all_targeted_comparisons.done"))


rule targeted_quast_comparison:
    input:
        # The input here is just a placeholder to connect to the checkpoint.
        # The actual files are found dynamically in the 'run' block.
        # This input ensures this rule only runs after the binner for the sample has run.
        binned_dirs=lambda wildcards: expand(
            os.path.join(ANNOTATION_DIR, wildcards.sample, "binned_contigs", "{assembler}"),
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        report=os.path.join(STATS_DIR, "targeted_quast", "{sample}", "{virus}", "report.html"),
        out_dir=directory(os.path.join(STATS_DIR, "targeted_quast", "{sample}", "{virus}")),
        transposed=os.path.join(STATS_DIR, "targeted_quast", "{sample}", "{virus}", "transposed_report.tsv")
    params:
        reference=lambda wildcards: config["viruses_of_interest"][wildcards.virus],
    log:
        os.path.join(LOG_DIR, "targeted_quast_comparison", "{sample}_{virus}.log")
    run:
        # At execution time, find all existing FASTA files for this sample/virus combo.
        fastas_to_compare = glob.glob(os.path.join(ANNOTATION_DIR, wildcards.sample, "binned_contigs", "*", f"{wildcards.virus}.fasta"))

        if fastas_to_compare:
            labels = sorted([path.split(os.sep)[-2] for path in fastas_to_compare])
            fastas_to_compare.sort() # Ensure labels and files are in the same order
            
            labels_str = ",".join(labels)
            fastas_str = " ".join(fastas_to_compare)

            shell(
                "quast.py -o {output.out_dir} -r {params.reference} "
                "-l '{labels_str}' {fastas_str} &> {log}"
            )
        else:
            # This logic should not be reached if the aggregator works correctly,
            # but it is kept for safety.
            shell("mkdir -p {output.out_dir} && touch {output.report} && touch {output.transposed}")


# Step 9: Calculate assembly statistics
rule calculate_stats:
    input:
        get_assembly_fasta
    output:
        os.path.join(STATS_DIR, "contig_stats", "{sample}_{assembler}.stats.txt")
    log:
        os.path.join(LOG_DIR, "calculate_stats", "{sample}_{assembler}.log")
    shell:
        "stats.sh {input} format=5 > {output} 2> {log}"

# Step 10: Map all QC'd reads back to each assembly
rule map_reads:
    input:
        contigs=get_assembly_fasta,
        reads=os.path.join(QC_DIR, "{sample}.qc.fastq")
    output:
        os.path.join(STATS_DIR, "reads_to_contigs", "{sample}_{assembler}.bam")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "map_reads", "{sample}_{assembler}.log")
    shell:
        "minimap2 -aY -t {threads} -x map-ont {input.contigs} {input.reads} 2> {log} | "
        "samtools sort -@ {threads} --output-fmt BAM -o {output}"

rule aggregate_contig_mappings:
    input:
        get_possible_contig_mappings
    output:
        touch(os.path.join(STATS_DIR, "all_contig_mappings.done"))

# Unidentified step: Map all contigs back to reference
rule map_binned_contigs_to_reference:
    input:
        binned_fasta=os.path.join(ANNOTATION_DIR, "{sample}", "binned_contigs", "{assembler}", "{virus}.fasta")
    output:
        bam=os.path.join(STATS_DIR, "contigs_to_ref", "{sample}", "{assembler}", "{virus}.bam")
    params:
        # The reference is looked up dynamically from the config using the {virus} wildcard.
        reference=lambda wildcards: config["viruses_of_interest"][wildcards.virus]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "map_binned_contigs", "{sample}_{assembler}_{virus}.log")
    shell:
        """
        minimap2 -ax asm5 -t {threads} {params.reference} {input.binned_fasta} \
        2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        """

# Step 11: Run CheckV on each assembly
rule run_checkv:
    input:
        get_assembly_fasta
    output:
        folder=directory(os.path.join(STATS_DIR, "checkv", "{sample}_{assembler}")),
        quality_summary=os.path.join(STATS_DIR, "checkv", "{sample}_{assembler}", "quality_summary.tsv")
    params:
        db=config["paths"]["checkv_db"]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "run_checkv", "{sample}_{assembler}.log")
    shell:
        "checkv end_to_end {input} {output.folder} -d {params.db} -t {threads} &> {log}"

# RULES FOR SUMMARY REPORT
rule summarize_sample_checkv:
    input:
        lambda wildcards: expand(
            os.path.join(STATS_DIR, "checkv", f"{wildcards.sample}_{{assembler}}", "quality_summary.tsv"),
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        checkv_csv=os.path.join(STATS_DIR, "per_sample", "{sample}_checkv_summary.csv"),
        miuvig_csv=os.path.join(STATS_DIR, "per_sample", "{sample}_miuvig_summary.csv")
    params:
        script="scripts/summarize_checkv.py"
    log:
        os.path.join(LOG_DIR, "summarize_checkv", "{sample}.log")
    shell:
        "python {params.script} {output.checkv_csv} {output.miuvig_csv} {input} > {log} 2>&1"

rule summarize_benchmarks:
    input:
        benchmarks=expand("benchmarks/{process}/{sample}.log", process=["classify_reads_diamond"] + [f"assemble_{asm}" for asm in ACTIVE_ASSEMBLERS], sample=SAMPLES),
        assembly_stats=expand(os.path.join(STATS_DIR, "contig_stats", "{sample}_{assembler}.stats.txt"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        bams=expand(os.path.join(STATS_DIR, "reads_to_contigs", "{sample}_{assembler}.bam"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS)
    output:
        benchmark_csv=os.path.join(RESULTS_DIR, "report", "benchmark_summary.csv"),
        assembly_csv=os.path.join(RESULTS_DIR, "report", "assembly_summary.csv"),
        per_sample_benchmarks=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_benchmark_summary.csv"), sample=SAMPLES),
        per_sample_assemblies=expand(os.path.join(STATS_DIR, "per_sample", "{sample}_assembly_summary.csv"), sample=SAMPLES)
    params:
        script="scripts/summarize_benchmarks.py",
        per_sample_dir=os.path.join(STATS_DIR, "per_sample")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "summarize_benchmarks.log")
    run:
        all_inputs = " ".join(input.benchmarks + input.assembly_stats + input.bams)
        
        shell(
            "python {params.script} --threads {threads} {all_inputs} > {log} 2>&1"
        )
        
        shell("mv benchmark_summary.csv {output.benchmark_csv}")
        shell("mv assembly_summary.csv {output.assembly_csv}")
        # Move per-sample files, if they exist, to the per_sample stats directory
        shell("mv *_benchmark_summary.csv {params.per_sample_dir}/ 2>/dev/null || true")
        shell("mv *_assembly_summary.csv {params.per_sample_dir}/ 2>/dev/null || true")
		
rule gather_versions:
     output:
        versions=os.path.join(RESULTS_DIR, "report", "versions.tsv")
     shell:
        """
        {{
        set -euo pipefail
        # Create a header
        echo "Name Version Channel"
        conda list --fields name,version,channel_name | grep -E "flye|canu|raven|plass|myloasm|metamdbg|wtdbg|shasta|miniasm|diamond" | awk '{{print $1, $2, $3}}'
        }} > {output.versions}
        """

# --- QUARTO BOOK GENERATION ---
include: "rules/generate_report.smk"
