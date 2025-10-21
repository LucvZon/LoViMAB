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
            os.path.join(LOG_DIR, "assemble_metaflye", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_metaflye", "{sample}.log")
        shell:
            "flye --nano-raw {input} --meta -o {output.dir} --threads {threads} &> {log}"

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
            os.path.join(LOG_DIR, "assemble_penguin", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_penguin", "{sample}.log")
        shell:
            "penguin nuclassemble {input} {output.fasta} {output.tmp_dir} "
            "--min-contig-len {params.min_len} --min-seq-id {params.min_id} "
            "--threads {threads} &> {log}"

if ASSEMBLERS_CONFIG.get("raven", False):
    rule assemble_raven:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            os.path.join(ASSEMBLY_DIR, "{sample}", "raven", "assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "assemble_raven", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_raven", "{sample}.log")
        shell:
            """
            raven --threads {threads} -p 2 {input} > {output} 2> {log}
            rm raven.cereal
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
            os.path.join(LOG_DIR, "assemble_canu", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_canu", "{sample}.log")
        shell:
            """
            canu -p canu_assembly \
            -d {output.dir} \
            genomeSize={params.genome_size} \
            maxThreads={threads} \
            maxMemory={params.memory} \
            -nanopore {input} \
            useGrid=false 2> {log}
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
            os.path.join(LOG_DIR, "assemble_myloasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_myloasm", "{sample}.log")
        shell:
            """
            myloasm {input} \
            -o {output.dir} \
            --min-reads-contig {params.min_reads} \
            --min-ol {params.min_overlap} \
            -t {threads} 2> {log}
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
            os.path.join(LOG_DIR, "assemble_metamdbg", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_metamdbg", "{sample}.log")
        shell:
            """
            metaMDBG asm --in-ont {input} \
            --out-dir {output.dir} \
            --min-read-overlap {params.min_overlap} \
            --min-read-identity {params.min_id} \
            --threads {threads} 2> {log}

            gzip --decompress -c {output.dir}/contigs.fasta.gz > {output.fasta}
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
            os.path.join(LOG_DIR, "assemble_wtdbg2", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_wtdbg2", "{sample}.log")
        shell:
            """
            mkdir -p {output.dir}

            wtdbg2 \
            -i results/2_read_classification/Sample1.target_reads.fastq \
            -o {output.dir}/dbg \
            -t {threads} \
            -x ont \
            -L {params.min_read_length} \
            -e 2 \
            --ctg-min-length {params.min_contig_length} 2> {log}
            
            wtpoa-cns -t {threads} -i {output.dir}/dbg.ctg.lay.gz -fo {output.fasta}
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
            os.path.join(LOG_DIR, "assemble_shasta", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_shasta", "{sample}.log")
        shell:
            """
            shasta \
            --input {input} \
            --threads {threads} \
            --config Nanopore-R10-Fast-Nov2022 \
            --Reads.minReadLength {params.min_read_length} \
            --assemblyDirectory {output.dir} \
            --Align.minAlignedMarkerCount 30 \
            --MarkerGraph.minCoverage 3 2> {log}
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
            os.path.join(LOG_DIR, "assemble_miniasm", "{sample}.log")
        benchmark:
            os.path.join(BENCH_DIR, "assemble_miniasm", "{sample}.log")
        shell:
            """
            # Create draft assembly graph
            minimap2 -t {threads} -x ava-ont {input} {input} > {output.dir}/overlaps.paf
            miniasm -s {params.min_overlap} -f {input} {output.dir}/overlaps.paf > {output.dir}/raw_assembly.gfa 2> {log}

            # Convert graph to fasta
            gfatools gfa2fa {output.dir}/raw_assembly.gfa > {output.dir}/raw_assembly.fasta

            # First polishing round
            minimap2 -t {threads} -x map-ont {output.dir}/raw_assembly.fasta {input} > {output.dir}/polished_overlaps_1.paf
            racon -t {threads} {input} {output.dir}/polished_overlaps_1.paf {output.dir}/raw_assembly.fasta > {output.dir}/polished_assembly_1.fasta
            # Second polishing round
            minimap2 -t {threads} -x map-ont {output.dir}/polished_assembly_1.fasta {input} > {output.dir}/polished_overlaps_2.paf
            racon -t {threads} {input} {output.dir}/polished_overlaps_2.paf {output.dir}/polished_assembly_1.fasta > {output.fasta}
            """
