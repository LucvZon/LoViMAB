# File: workflow/rules/assemblers.smk
# This file contains all rules related to genomic assembly.
# It expects variables like ASSEMBLERS_CONFIG, READ_CLASSIFICATION_DIR,
# ASSEMBLY_DIR, LOG_DIR, and BENCH_DIR to be defined in the main Snakefile.

# --- ASSEMBLY RULES ---

if ASSEMBLERS_CONFIG.get("metaflye", False):
    rule assemble_metaflye:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "metaflye")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "metaflye", "assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "metaflye", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "metaflye", "{sample}.log")
        shell:
            """
            # Wrap the flye command in parentheses and use '||' for a fallback.
            # The '&>' redirects both stdout and stderr to the log file.
            (flye --nano-raw {input} --meta --min-overlap 1000 -o {output.dir} --threads {threads} &> {log}) \
            || \
            (echo "Flye failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if ASSEMBLERS_CONFIG.get("penguin", False):
    rule assemble_penguin:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "penguin", "contigs.fasta"),
            tmp_dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "penguin", "temp_files"))
        params:
            min_len=config["params"]["penguin_min_contig_len"],
            min_id=config["params"]["penguin_min_seq_id"]
        threads:
            config["params"]["threads"]		
        log:
            os.path.join(LOG_DIR, "primary", "penguin", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "penguin", "{sample}.log")
        shell:
            """
            (penguin nuclassemble {input} {output.fasta} {output.tmp_dir} \
            --min-contig-len {params.min_len} --min-seq-id {params.min_id} \
            --threads {threads} &> {log}) \
            || \
            (echo "PenguiN failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             touch {output.fasta})
            """

if ASSEMBLERS_CONFIG.get("raven", False):
    rule assemble_raven:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            os.path.join(ASSEMBLY_DIR, "{sample}", "raven", "assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "raven", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "raven", "{sample}.log")
        shell:
            """
            (raven --threads {threads} -p 2 {input} > {output} 2> {log}
            rm raven.cereal) \
            || \
            (echo "Raven failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output})
            """

if ASSEMBLERS_CONFIG.get("canu", False):
    rule assemble_canu:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "canu")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "canu", "canu_assembly.contigs.fasta")
        params:
            genome_size=config["params"]["canu_genome_size"],
            # A good rule of thumb: 4GB of memory per thread for Canu
            memory=lambda wildcards, threads: threads * 4
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "canu", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "canu", "{sample}.log")
        shell:
            """
            # Wrap the canu command in parentheses and use '||' for a fallback
            (canu -p canu_assembly \
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

if ASSEMBLERS_CONFIG.get("myloasm", False):
    rule assemble_myloasm:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "myloasm")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "myloasm", "assembly_primary.fa")
        params:
            min_reads=config["params"]["myloasm_min_reads"],
            min_overlap=config["params"]["myloasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "myloasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "myloasm", "{sample}.log")
        shell:
            """
            # Try to run myloasm. If it fails (exits non-zero), run the command after '||'
            (myloasm {input} \
                -o {output.dir} \
                --min-reads-contig {params.min_reads} \
                --min-ol {params.min_overlap} \
                -t {threads} 2> {log}) \
            || \
            (echo "Myloasm failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            """

if ASSEMBLERS_CONFIG.get("metamdbg", False):
    rule assemble_metamdbg:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "metamdbg")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "metamdbg", "contigs.fasta")
        params:
            min_overlap=config["params"]["metamdbg_min_overlap"],
            min_id=config["params"]["metamdbg_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "metamdbg", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "metamdbg", "{sample}.log")
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

if ASSEMBLERS_CONFIG.get("wtdbg2", False):
    rule assemble_wtdbg2:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "wtdbg2")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "wtdbg2", "contigs.fasta")
        params:
            min_read_length=config["params"]["wtdbg2_min_read_len"],
            min_contig_length=config["params"]["wtdbg2_min_contig_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "wtdbg2", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "wtdbg2", "{sample}.log")
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

if ASSEMBLERS_CONFIG.get("shasta", False):
    rule assemble_shasta:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "shasta")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "shasta", "Assembly.fasta")
        params:
            min_read_length=config["params"]["shasta_min_read_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "shasta", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "shasta", "{sample}.log")
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

if ASSEMBLERS_CONFIG.get("miniasm", False):
    rule assemble_miniasm:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "miniasm")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "miniasm", "final_assembly.fasta")
        params:
            min_overlap=config["params"]["miniasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "miniasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "primary", "miniasm", "{sample}.log")
        shell:
            """
            # Create draft assembly graph
            minimap2 -t {threads} -x ava-ont {input} {input} > {output.dir}/overlaps.paf
            # Redirect miniasm stdout and stderr to log
            miniasm -s {params.min_overlap} -f {input} {output.dir}/overlaps.paf > {output.dir}/raw_assembly.gfa 2> {log}

            # Convert graph to fasta
            gfatools gfa2fa {output.dir}/raw_assembly.gfa > {output.dir}/raw_assembly.fasta

            # --- ROBUSTNESS CHECK ---
            # Check if the raw assembly FASTA is empty. The '-s' flag checks if a file has a size greater than zero.
            if [ -s {output.dir}/raw_assembly.fasta ]; then
                # If the file is NOT empty, proceed with polishing
                echo "Miniasm produced contigs. Proceeding with Racon polishing." >> {log}

                # First polishing round
                minimap2 -t {threads} -x map-ont {output.dir}/raw_assembly.fasta {input} > {output.dir}/polished_overlaps_1.paf 2>> {log}
                racon -t {threads} {input} {output.dir}/polished_overlaps_1.paf {output.dir}/raw_assembly.fasta > {output.dir}/polished_assembly_1.fasta 2>> {log}
                
                # Second polishing round
                minimap2 -t {threads} -x map-ont {output.dir}/polished_assembly_1.fasta {input} > {output.dir}/polished_overlaps_2.paf 2>> {log}
                racon -t {threads} {input} {output.dir}/polished_overlaps_2.paf {output.dir}/polished_assembly_1.fasta > {output.fasta} 2>> {log}

            else
                # If the file IS empty, skip polishing and create an empty final file
                echo "Miniasm produced no contigs for sample {wildcards.sample}. Skipping polishing." >> {log}
                touch {output.fasta}
            fi
            """
			
if ASSEMBLERS_CONFIG.get("cap3", False):
    rule assemble_cap3:
        input:
            # Use expand to gather contigs from ALL other active assemblers for this sample.
            # PRIMARY_ASSEMBLERS is the list defined in the main snakefile.
            expand(os.path.join(ASSEMBLY_DIR, "{{sample}}", "{assembler}", "renamed_contigs.fasta"), assembler=PRIMARY_ASSEMBLERS)
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "cap3")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "cap3", "final_reassembly.fasta")
        params:
            overlap=config["params"]["cap3_overlap_len"],
            identity=config["params"]["cap3_identity_pct"],
            combined_fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "cap3", "combined_contigs.fasta")
        log:
            os.path.join(LOG_DIR, "assemble_cap3", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_cap3", "{sample}.log")
        shell:
            """
            # Combine all input contigs from the primary assemblers into a single file
            cat {input} > {params.combined_fasta}

            # Run cap3 on the combined contigs. It outputs files based on the input filename.
            # We redirect all screen output to the log file.
            cap3 {params.combined_fasta} -o {params.overlap} -p {params.identity} > {log} 2>&1

            # CAP3 creates a file named '{params.combined_fasta}.cap.contigs'.
            # We move and rename it to our desired final output name.
            mv {params.combined_fasta}.cap.contigs {output.fasta}
            """
