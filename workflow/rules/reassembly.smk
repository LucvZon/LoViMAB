# This snakefile performs reassembly with each active assembler

# --- Helper function to get correct assembly output file ---
def get_reassembly_fasta(wildcards):
    assembler = wildcards.assembler
    if assembler == "metaflye":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "penguin":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "raven":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "canu":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "canu_assembly.contigs.fasta")
    elif assembler == "myloasm":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly_primary.fa")
    elif assembler == "metamdbg":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "wtdbg2":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "shasta":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "Assembly.fasta")
    elif assembler == "miniasm":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "final_assembly.fasta")
    elif assembler == "cap3":
         return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "final_reassembly.fasta")
    # Add other assemblers here if needed in the future

# Place assembly rules here
if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_metaflye:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "metaflye")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "metaflye", "assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_metaflye", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_metaflye", "{sample}.log")
        shell:
            """
            (flye --nano-raw {input} --meta --min-overlap 1000 -o {output.dir} --threads {threads} &> {log}) \
            || \
            (echo "Flye failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_penguin:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "penguin", "contigs.fasta"),
            tmp_dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "penguin", "temp_files"))
        params:
            min_len=config["params"]["penguin_min_contig_len"],
            min_id=config["params"]["penguin_min_seq_id"]
        threads:
            config["params"]["threads"]		
        log:
            os.path.join(LOG_DIR, "reassemble_penguin", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_penguin", "{sample}.log")
        shell:
            """
            (penguin nuclassemble {input} {output.fasta} {output.tmp_dir} \
            --min-contig-len {params.min_len} --min-seq-id {params.min_id} \
            --threads {threads} &> {log}) \
            || \
            (echo "PenguiN failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_raven:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            os.path.join(REASSEMBLY_DIR, "{sample}", "raven", "assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_raven", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_raven", "{sample}.log")
        shell:
            """
            (raven --threads {threads} -p 2 {input} > {output} 2> {log}
            rm raven.cereal) \
            || \
            (echo "Raven failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_canu:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "canu")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "canu", "canu_assembly.contigs.fasta")
        params:
            genome_size=config["params"]["canu_genome_size"],
            # A good rule of thumb: 4GB of memory per thread for Canu
            memory=lambda wildcards, threads: threads * 4
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_canu", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_canu", "{sample}.log")
        shell:
            """
            (canu -assemble -corrected \
            -p canu_assembly \
            -d {output.dir} \
            genomeSize={params.genome_size} \
            maxThreads={threads} \
            maxMemory={params.memory} \
            -nanopore {input} \
            useGrid=false 2> {log}) \
            || \
            (echo "Canu failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_myloasm:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "myloasm")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "myloasm", "assembly_primary.fa")
        params:
            min_reads=config["params"]["myloasm_min_reads"],
            min_overlap=config["params"]["myloasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_myloasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_myloasm", "{sample}.log")
        shell:
            """
            (myloasm {input} \
            -o {output.dir} \
            --min-reads-contig 1 \
            --min-ol {params.min_overlap} \
            -t {threads} 2> {log}) \
            || \
            (echo "Myloasm failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_metamdbg:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "metamdbg")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "metamdbg", "contigs.fasta")
        params:
            min_overlap=config["params"]["metamdbg_min_overlap"],
            min_id=config["params"]["metamdbg_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_metamdbg", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_metamdbg", "{sample}.log")
        shell:
            """
            (metaMDBG asm --in-ont {input} \
            --out-dir {output.dir} \
            --min-read-overlap {params.min_overlap} \
            --min-read-identity {params.min_id} \
            --threads {threads} 2> {log}

            gzip --decompress -c {output.dir}/contigs.fasta.gz > {output.fasta}) \
            || \
            (echo "metaMDBG failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_wtdbg2:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "wtdbg2")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "wtdbg2", "contigs.fasta")
        params:
            min_read_length=config["params"]["wtdbg2_min_read_len"],
            min_contig_length=config["params"]["wtdbg2_min_contig_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_wtdbg2", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_wtdbg2", "{sample}.log")
        shell:
            """
            mkdir -p {output.dir}

            (wtdbg2 \
            -i {input} \
            -o {output.dir}/dbg \
            -t {threads} \
            -x ont \
            -L {params.min_read_length} \
            -e 2 \
            --ctg-min-length {params.min_contig_length} 2> {log}
            
            wtpoa-cns -t {threads} -i {output.dir}/dbg.ctg.lay.gz -fo {output.fasta}) \
            || \
            (echo "Wtdbg2 failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_shasta:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "shasta")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "shasta", "Assembly.fasta")
        params:
            min_read_length=config["params"]["shasta_min_read_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_shasta", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_shasta", "{sample}.log")
        shell:
            """
            (shasta \
            --input {input} \
            --threads {threads} \
            --config Nanopore-R10-Fast-Nov2022 \
            --Reads.minReadLength {params.min_read_length} \
            --assemblyDirectory {output.dir} \
            --Align.minAlignedMarkerCount 30 \
            --MarkerGraph.minCoverage 3 2> {log}) \
            || \
            (echo "Shasta failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_miniasm:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "miniasm")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "miniasm", "final_assembly.fasta")
        params:
            min_overlap=config["params"]["miniasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "reassemble_miniasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "reassemble_miniasm", "{sample}.log")
        shell:
            """
            # Create draft assembly graph
            minimap2 -t {threads} -x ava-ont {input} {input} > {output.dir}/overlaps.paf
            miniasm -s {params.min_overlap} -f {input} {output.dir}/overlaps.paf > {output.dir}/raw_assembly.gfa 2> {log}

            # Convert graph to fasta
            gfatools gfa2fa {output.dir}/raw_assembly.gfa > {output.fasta}
            """


# Place quality assessment rules here
# Step 11: Calculate assembly statistics
rule reassembly_calculate_stats:
    input:
        get_reassembly_fasta
    output:
        os.path.join("results", "6_reassembly_stats", "contig_stats", "{sample}_{assembler}.stats.txt")
    log:
        os.path.join(LOG_DIR, "reassembly_calculate_stats", "{sample}_{assembler}.log")
    shell:
        "stats.sh {input} format=5 > {output} 2> {log}"

# Map QC'd reads back to reassembled contigs
rule map_reads_to_reassembly:
    input:
        contigs=get_reassembly_fasta,
        reads=os.path.join(QC_DIR, "{sample}.qc.fastq")
    output:
        os.path.join("results", "6_reassembly_stats", "reads_to_contigs", "{sample}_{assembler}.bam")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "reassemble_map_reads", "{sample}_{assembler}.log")
    shell:
        "minimap2 -aY -t {threads} -x map-ont {input.contigs} {input.reads} 2> {log} | "
        "samtools sort -@ {threads} --output-fmt BAM -o {output}"

# Map original contigs to reassembled contigs
rule map_contigs_to_reassembly:
    input:
        reassembled=get_reassembly_fasta,
        original=get_assembly_fasta
    output:
        os.path.join("results", "6_reassembly_stats", "original_to_reassembled", "{sample}_{assembler}.bam")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "contigs_to_reassembly", "{sample}_{assembler}.log")
    shell:
        "minimap2 -aY -t {threads} -x map-ont {input.reassembled} {input.original} 2> {log} | "
        "samtools sort -@ {threads} --output-fmt BAM -o {output}"

# Utilize old summary scripts logic
rule summarize_reassembly_benchmarks:
    input:
        benchmarks=expand("benchmarks/{process}/{sample}.log", process=[f"reassemble_{asm}" for asm in ACTIVE_ASSEMBLERS], sample=SAMPLES),
        assembly_stats=expand(os.path.join("results", "6_reassembly_stats", "contig_stats", "{sample}_{assembler}.stats.txt"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS),
        bams=expand(os.path.join("results", "6_reassembly_stats", "reads_to_contigs", "{sample}_{assembler}.bam"), sample=SAMPLES, assembler=ACTIVE_ASSEMBLERS)
    output:
        benchmark_csv=os.path.join(RESULTS_DIR, "report", "reassembly_benchmark_summary.csv"),
        assembly_csv=os.path.join(RESULTS_DIR, "report", "reassembly_assembly_summary.csv"),
        per_sample_benchmarks=expand(os.path.join("results", "6_reassembly_stats", "per_sample", "{sample}_benchmark_summary.csv"), sample=SAMPLES),
        per_sample_assemblies=expand(os.path.join("results", "6_reassembly_stats", "per_sample", "{sample}_assembly_summary.csv"), sample=SAMPLES)
    params:
        script="scripts/summarize_benchmarks.py",
        per_sample_dir=os.path.join("results", "6_reassembly_stats", "per_sample")
    threads:
        config["params"]["threads"]
    log:
        os.path.join(LOG_DIR, "reassembly_summarize_benchmarks.log")
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

