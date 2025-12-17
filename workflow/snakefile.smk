# ===================================================================
#          SNAKEMAKE WORKFLOW FOR TARGETED METAGENOMIC ASSEMBLY
# ===================================================================

# VERSION: 1.0.0

import os
import glob

# --- Load Configuration ---
configfile: "config/lovimab.yaml"

# --- Global Variables ---
SAMPLES = list(config["samples"].keys())
ASSEMBLERS_CONFIG = config["assemblers"]
REASSEMBLY_CONFIG = config["reassembly"]
ACTIVE_ASSEMBLERS = [asm for asm, active in ASSEMBLERS_CONFIG.items() if active]
PRIMARY_ASSEMBLERS = [asm for asm in ACTIVE_ASSEMBLERS if asm != "cap3"]
ASSEMBLY_TYPES = ["primary", "secondary", "final"]

# --- Define output directories ---
RESULTS_DIR = "results"
QC_DIR = os.path.join(RESULTS_DIR, "1_quality_control")
READ_CLASSIFICATION_DIR = os.path.join(RESULTS_DIR, "2_read_classification")
ASSEMBLY_DIR = os.path.join(RESULTS_DIR, "3_assemblies")
REASSEMBLY_DIR = os.path.join(RESULTS_DIR, "4_reassemblies")
CLUSTER_DIR = os.path.join(RESULTS_DIR, "5_clusters")
FINAL_DIR = os.path.join(RESULTS_DIR, "6_final_assembly")
ANNOTATION_DIR = os.path.join(RESULTS_DIR, "7_annotation")
STATS_DIR = os.path.join(RESULTS_DIR, "8_stats_and_qc")
LOG_DIR = "logs"
BENCH_DIR = "benchmarks"

# --- Handle Optional Viruses of Interest ---
# Safely get the dictionary, defaulting to empty if missing or None
VIRUSES_OF_INTEREST = config.get("viruses_of_interest")
if VIRUSES_OF_INTEREST is None:
    VIRUSES_OF_INTEREST = {}

# --- Onstart: Validate Configuration ---
onstart:
    # Check if every file path listed in configfile exists.
    # We now use the safe VIRUSES_OF_INTEREST variable
    for virus_name, ref_path in VIRUSES_OF_INTEREST.items():
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
    elif assembler == "cap3":
         return os.path.join(ASSEMBLY_DIR, wildcards.sample, assembler, "final_reassembly.fasta")
    # Add other assemblers here if needed in the future

# --- Helper function to handle different assembly types ---
def get_assembly_by_type(wildcards):
    """
    Returns the path to the FASTA file based on the assembly_type wildcard.
    """
    a_type = wildcards.assembly_type
    sample = wildcards.sample
    assembler = wildcards.assembler

    if a_type == "primary":
        # Uses your existing logic for primary assemblies
        return get_assembly_fasta(wildcards)
        
    elif a_type == "secondary":
        # Uses your existing logic for secondary (reassembly)
        # Ensure this function is available here or import/define it
        return get_reassembly_fasta(wildcards)
        
    elif a_type == "final":
        # The new path for the mmseqs-merged assembly
        return os.path.join(FINAL_DIR, sample, assembler, "final_assembly.fasta")
        
    return "UNKNOWN_ASSEMBLY_TYPE"

# ===================================================================
#                              TARGET RULE
# ===================================================================
rule all:
    input:
        # --- Generic results ---
        expand(os.path.join(STATS_DIR, "contig_stats", "{assembly_type}", "{sample}_{assembler}.stats.txt"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type=ASSEMBLY_TYPES),
        expand(os.path.join(STATS_DIR, "reads_to_contigs", "{assembly_type}", "{sample}_{assembler}.bam"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type=ASSEMBLY_TYPES),
        expand(os.path.join(STATS_DIR, "checkv", "{assembly_type}", "{sample}_{assembler}"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type=ASSEMBLY_TYPES),
        expand(os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type=ASSEMBLY_TYPES),
        
        # Temporary rule to trigger rule combine_assemblies
        # expand(os.path.join("results", "6_reassembly_stats", "contig_stats", "{sample}_{assembler}.stats.txt"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),

        # Some temporary outputs until we fix the final report, then these will become obsolete
        # BENCHMARK/ASSEMBLY-STATS
        expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "benchmark_summary.csv"), assembly_type=["primary", "secondary", "final"]),
        expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_benchmark_summary.csv"), assembly_type=["primary", "secondary", "final"], sample=SAMPLES),
        # CHECKV
        expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_checkv_summary.csv"), assembly_type=["primary", "secondary", "final"], sample=SAMPLES),
        expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_miuvig_summary.csv"), assembly_type=["primary", "secondary", "final"], sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "checkv_global_summary.csv"), assembly_type=["primary", "secondary", "final"]),
        expand(os.path.join(RESULTS_DIR, "report", "{assembly_type}", "miuvig_global_summary.csv"), assembly_type=["primary", "secondary", "final"]),

        # --- Specific, Targeted Reports ---
        os.path.join(STATS_DIR, "all_targeted_comparisons.done"),
        os.path.join(STATS_DIR, "all_contig_mappings.done"),
        os.path.join(RESULTS_DIR, "report", "final_summary_report", "index.html")


# ===================================================================
#                          WORKFLOW RULES
# ===================================================================

def get_possible_quast_comparisons(wildcards):
    """
    Scans checkpoint outputs for all assembly types, samples, and assemblers
    to determine which targeted QUAST reports to generate.
    """
    # Return empty immediately if no viruses of interest are configured
    if not VIRUSES_OF_INTEREST:
        return []

    found_combos = set()

    # Iterate through every combination    
    for a_type in ASSEMBLY_TYPES:
        for sample in SAMPLES:
            for assembler in ACTIVE_ASSEMBLERS:
                # Get the specific checkpoint output for this Type/Sample/Assembler combo
                binned_dir = checkpoints.bin_contigs_by_taxonomy.get(
                    assembly_type=a_type, sample=sample, assembler=assembler
                ).output.binned_dir
                
                # Scan for fasta files
                # Path: results/4_annotation/{type}/{sample}/binned_contigs/{assembler}/*.fasta
                for fasta_file in glob.glob(os.path.join(str(binned_dir), "*.fasta")):
                    virus_name = os.path.basename(fasta_file).replace(".fasta", "")
                    
                    if virus_name in VIRUSES_OF_INTEREST:
                        # Store the unique (type, sample, virus) combination
                        found_combos.add((a_type, sample, virus_name))

    # Build final paths
    possible_reports = []
    for a_type, sample, virus_name in found_combos:
        possible_reports.append(
            os.path.join(STATS_DIR, "targeted_quast", a_type, sample, virus_name, "report.html")
        )
        
    return possible_reports

def get_possible_contig_mappings(wildcards):
    """
    This function scans the checkpoint output to determine which binned contigs
    can be mapped to their corresponding reference of interest.
    """
    # Return empty immediately if no viruses of interest are configured
    if not VIRUSES_OF_INTEREST:
        return []

    possible_bam_files = []
    
    for a_type in ASSEMBLY_TYPES:
        for sample in SAMPLES:
            for assembler in ACTIVE_ASSEMBLERS:
                binned_dir = checkpoints.bin_contigs_by_taxonomy.get(
                    assembly_type=a_type, sample=sample, assembler=assembler
                ).output.binned_dir
                
                for fasta_file in glob.glob(os.path.join(str(binned_dir), "*.fasta")):
                    virus_name = os.path.basename(fasta_file).replace(".fasta", "")
                    
                    if virus_name in VIRUSES_OF_INTEREST:
                        possible_bam_files.append(
                            os.path.join(STATS_DIR, "contigs_to_ref", a_type, sample, assembler, f"{virus_name}.bam")
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

# Step 2: Trim adapter sequences 
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

# Step 3: Perform quality control on merged reads
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

# Step 4: Classify reads with DIAMOND against a custom database
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
        os.path.join(LOG_DIR, "primary", "classify_reads_diamond", "{sample}.log")
    benchmark:
        os.path.join(BENCH_DIR, "primary", "classify_reads_diamond", "{sample}.log")
    shell:
        "diamond blastx {params.sensitivity_flag} -d {params.db} -q {input} -o {output} "
        "-f 6 qseqid -k 1 --threads {threads} &> {log}" # Using -k 1 for best hit

# Step 5: Extract target reads based on DIAMOND classification
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
		
# Step 5b: Calculate the percentage of reads kept after classification
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

# Step 5c: calculate average read length
rule calculate_mean_length:
    input:
        target_reads_file=os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
    output:
        os.path.join(STATS_DIR, "per_sample", "{sample}_mean_read_length.txt")
    log:
        os.path.join(LOG_DIR, "calculate_mean_readlength", "{sample}.log")
    shell:
        """
        awk '{{if(NR%4==2){{count++; bases += length}}}} END{{if (count > 0) printf "%.0f", bases/count; else print "0"}}' {input.target_reads_file} > {output} 2> {log}
        """

# --- ASSEMBLY RULES ---
include: "rules/assemblers.smk"

# run this every time for now
rule combine_assemblies:
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY_DIR, wildcards.sample, "{assembler}", "renamed_contigs.fasta"), assembler=ACTIVE_ASSEMBLERS)
    output:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
    log:
        os.path.join(LOG_DIR, "combine_assemblies", "{sample}.log")
    shell:
        """
        # Combine all contigs from active assemblers into a single file per sample.
        cat {input} > {output.fasta} 2> {log}
        """

# --- REASSEMBLY RULES ---
include: "rules/reassembly.smk"

# --- POST-ASSEMBLY PROCESSING ---

# Step 6a: Rename contigs to ensure uniqueness before aggregation
rule rename_contigs_generic:
    input:
        fasta=get_assembly_by_type
    output:
        os.path.join(ANNOTATION_DIR, "staging", "{assembly_type}", "{sample}", "{assembler}.fasta")
    params:
        prefix="{assembler}_"
    log:
        os.path.join(LOG_DIR, "rename_contigs", "{assembly_type}", "{sample}_{assembler}.log")
    shell:
        "seqkit replace -p '^' -r '{params.prefix}' {input.fasta} > {output} 2> {log}"

# Step 6b: Aggregate all renamed contigs from a single sample
rule aggregate_contigs_generic:
    input:
        # Look into the staging directory for all active assemblers
        lambda wildcards: expand(
            os.path.join(ANNOTATION_DIR, "staging", "{assembly_type}", "{sample}", "{assembler}.fasta"),
            assembly_type=wildcards.assembly_type,
            sample=wildcards.sample,
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        # Output structure: results/4_annotation/{type}/{sample}/aggregated_contigs.fasta
        os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "aggregated_contigs.fasta")
    log:
        os.path.join(LOG_DIR, "aggregate_contigs", "{assembly_type}", "{sample}.log")
    shell:
        "cat {input} > {output} 2> {log}"

# Step 6c: Annotate the single aggregated contig file
rule annotate_aggregated_contigs_generic:
    input:
        os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "aggregated_contigs.fasta")
    output:
        os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "aggregated_annotation.tsv")
    params:
        db=config["paths"]["diamond_db"]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "annotate_aggregated", "{assembly_type}", "{sample}.log")
    shell:
        """
        # Check if the aggregated input FASTA is non-empty
        if [ -s {input} ]; then
            diamond blastx -d {params.db} -q {input} -o {output} \
                -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
                --threads {threads} -b 10 -c 1 &> {log}
        else
            echo "Aggregated contigs file is empty. Skipping DIAMOND annotation." > {log}
            touch {output}
        fi
        """

# Step 7: Split the aggregated annotation file back into per-assembler files
rule split_annotations_generic:
    input:
        os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "aggregated_annotation.tsv")
    output:
        os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "split", "{assembler}_annotation.tsv")
    params:
        prefix="{assembler}_"
    shell:
        "grep -E '^{params.prefix}' {input} > {output} || touch {output}"

# Step 8: Post-process the per-assembler annotation file
rule post_process_annotation_generic:
    input:
        annotation=os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "split", "{assembler}_annotation.tsv"),
        # Note: We use the STAGING file here as the source of contig sequences
        contigs=os.path.join(ANNOTATION_DIR, "staging", "{assembly_type}", "{sample}", "{assembler}.fasta")
    output:
        annotated=os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"),
        unannotated=os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "post_processed", "{assembler}_unannotated_contigs.tsv")
    params:
        script=workflow.source_path("../scripts/post_process_diamond_v1.0.py")
    log:
        os.path.join(LOG_DIR, "post_process", "{assembly_type}", "{sample}_{assembler}.log")
    shell:
        "python {params.script} -i {input.annotation} -c {input.contigs} "
        "-o {output.annotated} -u {output.unannotated} -log {log}"

# Step 9: Separate contigs based on annotation
checkpoint bin_contigs_by_taxonomy:
    input:
        annotation=os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "post_processed", "{assembler}_annotated_contigs.tsv"),
        contigs=os.path.join(ANNOTATION_DIR, "staging", "{assembly_type}", "{sample}", "{assembler}.fasta")
    output:
        binned_dir=directory(os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "binned_contigs", "{assembler}"))
    params:
        species_col=14, # Column name containing species information
        contig_col=1
    log:
        os.path.join(LOG_DIR, "bin_contigs", "{assembly_type}", "{sample}_{assembler}.log")
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

# Step 10: Run QUAST for all viruses of interest
rule targeted_quast_comparison:
    input:
        # The input here is just a placeholder to connect to the checkpoint.
        # The actual files are found dynamically in the 'run' block.
        # This input ensures this rule only runs after the binner for the sample has run.
        binned_dirs=lambda wildcards: expand(
            os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "binned_contigs", "{assembler}"),
            assembly_type=wildcards.assembly_type,
            sample=wildcards.sample,
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        report=os.path.join(STATS_DIR, "targeted_quast", "{assembly_type}", "{sample}", "{virus}", "report.html"),
        out_dir=directory(os.path.join(STATS_DIR, "targeted_quast", "{assembly_type}",  "{sample}", "{virus}")),
        transposed=os.path.join(STATS_DIR, "targeted_quast", "{assembly_type}", "{sample}", "{virus}", "transposed_report.tsv")
    params:
        reference=lambda wildcards: config["viruses_of_interest"][wildcards.virus]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "targeted_quast_comparison", "{assembly_type}", "{sample}_{virus}.log")
    run:
        # Find all binned fasta files for this specific Type/Sample/Virus combination
        # Path: results/4_annotation/{type}/{sample}/binned_contigs/*/{virus}.fasta
        search_pattern = os.path.join(
            ANNOTATION_DIR, 
            wildcards.assembly_type, 
            wildcards.sample, 
            "binned_contigs", 
            "*", 
            f"{wildcards.virus}.fasta"
        )
        fastas_to_compare = glob.glob(search_pattern)

        if fastas_to_compare:
            # Extract assembler names for labels (parent directory name)
            labels = sorted([path.split(os.sep)[-2] for path in fastas_to_compare])
            fastas_to_compare.sort() 
            
            labels_str = ",".join(labels)
            fastas_str = " ".join(fastas_to_compare)

            shell(
                "quast.py -o {output.out_dir} -r {params.reference} -t {threads} "
                "-l '{labels_str}' {fastas_str} &> {log}"
            )
        else:
            # Fallback
            shell("mkdir -p {output.out_dir} && touch {output.report} && touch {output.transposed}")


# Step 11: Calculate assembly statistics
rule calculate_stats_generic:
    input:
        # get_assembly_fasta (outdated :D)
        get_assembly_by_type
    output:
        os.path.join(STATS_DIR, "contig_stats", "{assembly_type}", "{sample}_{assembler}.stats.txt")
    log:
        os.path.join(LOG_DIR, "calculate_stats", "{assembly_type}", "{sample}_{assembler}.log")
    shell:
        "stats.sh {input} format=5 > {output} 2> {log}"

# Step 12: Map all QC'd reads back to each assembly
rule map_reads_generic:
    input:
        contigs=get_assembly_by_type,
        reads=os.path.join(QC_DIR, "{sample}.qc.fastq")
    output:
        os.path.join(STATS_DIR, "reads_to_contigs", "{assembly_type}", "{sample}_{assembler}.bam")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "map_reads", "{assembly_type}", "{sample}_{assembler}.log")
    shell:
        "minimap2 -aY -t {threads} -x map-ont {input.contigs} {input.reads} 2> {log} | "
        "samtools sort -@ {threads} --output-fmt BAM -o {output}"

rule aggregate_contig_mappings:
    input:
        get_possible_contig_mappings
    output:
        touch(os.path.join(STATS_DIR, "all_contig_mappings.done"))

# Step 13: Map binned contigs back to their virus of interest
rule map_binned_contigs_to_reference:
    input:
        binned_fasta=os.path.join(ANNOTATION_DIR, "{assembly_type}", "{sample}", "binned_contigs", "{assembler}", "{virus}.fasta")
    output:
        bam=os.path.join(STATS_DIR, "contigs_to_ref", "{assembly_type}", "{sample}", "{assembler}", "{virus}.bam")
    params:
        # The reference is looked up dynamically from the config using the {virus} wildcard.
        reference=lambda wildcards: config["viruses_of_interest"][wildcards.virus]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "map_binned_contigs", "{assembly_type}", "{sample}_{assembler}_{virus}.log")
    shell:
        """
        minimap2 -ax asm5 -t {threads} {params.reference} {input.binned_fasta} \
        2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        """

# Step 14: Run CheckV on each assembly
rule run_checkv_generic:
    input:
        get_assembly_by_type
    output:
        folder=directory(os.path.join(STATS_DIR, "checkv", "{assembly_type}", "{sample}_{assembler}")),
        quality_summary=os.path.join(STATS_DIR, "checkv", "{assembly_type}", "{sample}_{assembler}", "quality_summary.tsv")
    params:
        db=config["paths"]["checkv_db"]
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "run_checkv", "{assembly_type}", "{sample}_{assembler}.log")
    shell:
        """
        # First, check if the input FASTA file is non-empty
        if [ -s {input} ]; then
            # If it's not empty, run checkv normally
            checkv end_to_end {input} {output.folder} -d {params.db} -t {threads} &> {log}
        else
            # If it is empty, log the reason and create empty but valid output files
            echo "Input assembly for {wildcards.sample} / {wildcards.assembler} is empty. Skipping CheckV and creating dummy output files." > {log}
            
            # Create the main output directory
            mkdir -p {output.folder}
            
            # Create empty versions of all expected CheckV output files
            touch {output.folder}/proviruses.fna
            touch {output.folder}/viruses.fna
            
            # For TSV files, create them with just a header row to be safe
            echo -e "contig_id\\tcontig_length\\tviral_length\\taai_expected_length\\taai_completeness\\taai_confidence\\taai_error\\taai_num_hits\\taai_top_hit\\taai_id\\taai_af\\thmm_completeness_lower\\thmm_completeness_upper\\thmm_num_hits\\tkmer_freq" > {output.folder}/completeness.tsv
            echo -e "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings" > {output.folder}/quality_summary.tsv
            echo -e "contig_id\\tcontig_length\\ttotal_genes\\tviral_genes\\thost_genes\\tprovirus\\tproviral_length\\thost_length\\tregion_types\\tregion_lengths\\tregion_coords_bp\\tregion_coords_genes\\tregion_viral_genes\\tregion_host_genes" > {output.folder}/contamination.tsv
            echo -e "contig_id\\tcontig_length\\tkmer_freq\\tprediction_type\\tconfidence_level\\tconfidence_reason\\trepeat_length\\trepeat_count\\trepeat_n_freq\\trepeat_mode_base_freq\\trepeat_seq" > {output.folder}/complete_genomes.tsv
        fi
        """

# RULES FOR SUMMARY REPORT
rule summarize_sample_checkv:
    input:
        lambda wildcards: expand(
            os.path.join(STATS_DIR, "checkv", "{assembly_type}", "{sample}_{assembler}", "quality_summary.tsv"),
            assembly_type=wildcards.assembly_type,
            sample=wildcards.sample,
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        checkv_csv=os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_checkv_summary.csv"),
        miuvig_csv=os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_miuvig_summary.csv")
    params:
        script="scripts/summarize_checkv.py"
    log:
        os.path.join(LOG_DIR, "summarize_checkv", "{assembly_type}", "{sample}.log")
    shell:
        "python {params.script} {output.checkv_csv} {output.miuvig_csv} {input} > {log} 2>&1"

rule summarize_global_checkv:
    input:
        # Collect all quality summaries from all samples/assemblers
        lambda wildcards: expand(
            os.path.join(STATS_DIR, "checkv", "{assembly_type}", "{sample}_{assembler}", "quality_summary.tsv"),
            assembly_type=ASSEMBLY_TYPES,
            sample=SAMPLES,
            assembler=ACTIVE_ASSEMBLERS
        )
    output:
        checkv_csv=os.path.join(RESULTS_DIR, "report", "{assembly_type}", "checkv_global_summary.csv"),
        miuvig_csv=os.path.join(RESULTS_DIR, "report", "{assembly_type}", "miuvig_global_summary.csv")
    params:
        script="scripts/summarize_checkv.py"
    log:
        os.path.join(LOG_DIR, "summarize_checkv_global_{assembly_type}.log")
    shell:
        """
        # The script aggregates all input TSVs.
        # Since the input files now span ALL samples, the resulting CSVs will contain one row 
        # per assembler, summing the counts across all samples.
        python {params.script} {output.checkv_csv} {output.miuvig_csv} {input} > {log} 2>&1
        """

rule summarize_benchmarks_generic:
    input:
        # Use a lambda to dynamically pick the process list based on the wildcard
        benchmarks=lambda w: expand(
            os.path.join(BENCH_DIR, "{assembly_type}", "{process}", "{sample}.log"),
            assembly_type=w.assembly_type,
            sample=SAMPLES,
            # LOGIC: Add diamond only if type is 'primary', otherwise just use assemblers
            process=(["classify_reads_diamond"] + ACTIVE_ASSEMBLERS) if w.assembly_type == "primary" else ACTIVE_ASSEMBLERS
        ),
        # Fix: Pass wildcard string "{assembly_type}" to expand, do NOT use os.path.join keywords
        assembly_stats=expand(
            os.path.join(STATS_DIR, "contig_stats", "{assembly_type}", "{sample}_{assembler}.stats.txt"), 
            sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type="{assembly_type}"
        ),
        bams=expand(
            os.path.join(STATS_DIR, "reads_to_contigs", "{assembly_type}", "{sample}_{assembler}.bam"), 
            sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS, assembly_type="{assembly_type}"
        )
    output:
        # Fix: Removed invalid keywords inside os.path.join
        benchmark_csv=os.path.join(RESULTS_DIR, "report", "{assembly_type}", "benchmark_summary.csv"),
        assembly_csv=os.path.join(RESULTS_DIR, "report", "{assembly_type}", "assembly_summary.csv"),
        # Per-sample files
        per_sample_benchmarks=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_benchmark_summary.csv"), sample=SAMPLES, assembly_type="{assembly_type}"),
        per_sample_assemblies=expand(os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_assembly_summary.csv"), sample=SAMPLES, assembly_type="{assembly_type}")
    params:
        script="scripts/summarize_benchmarks.py",
        out_global_dir=os.path.join(RESULTS_DIR, "report", "{assembly_type}"),
        out_sample_dir=os.path.join(STATS_DIR, "per_sample", "{assembly_type}")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "summarize_benchmarks_{assembly_type}.log")
    run:
        import shutil
        import tempfile
        
        # 1. Create absolute paths for inputs/outputs/script because we will change CWD
        abs_script = os.path.abspath(params.script)
        abs_inputs = [os.path.abspath(f) for f in (input.benchmarks + input.assembly_stats + input.bams)]
        abs_out_global = os.path.abspath(params.out_global_dir)
        abs_out_sample = os.path.abspath(params.out_sample_dir)

        # 2. Create output directories if they don't exist
        os.makedirs(abs_out_global, exist_ok=True)
        os.makedirs(abs_out_sample, exist_ok=True)

        # 3. Create a temporary directory to run the script in
        # This prevents file collisions if 'primary' and 'secondary' run at the same time
        with tempfile.TemporaryDirectory() as temp_dir:
            
            # 4. Construct the command to run inside the temp dir
            # Note: We pass the absolute paths of inputs to the script
            cmd = f"cd {temp_dir} && python {abs_script} --threads {threads} {' '.join(abs_inputs)}"
            
            # 5. Run the command, logging to the original log file location
            shell(f"({cmd}) > {os.path.abspath(str(log))} 2>&1")
            
            # 6. Move the generated files from temp_dir to final destinations
            # Move Global Summaries
            shell(f"mv {temp_dir}/benchmark_summary.csv {output.benchmark_csv}")
            shell(f"mv {temp_dir}/assembly_summary.csv {output.assembly_csv}")
            
            # Move Per-Sample Summaries
            # We use '|| true' to avoid failure if the script didn't generate a specific sample file (e.g. empty inputs)
            shell(f"mv {temp_dir}/*_benchmark_summary.csv {abs_out_sample}/ 2>/dev/null || true")
            shell(f"mv {temp_dir}/*_assembly_summary.csv {abs_out_sample}/ 2>/dev/null || true")
		
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
