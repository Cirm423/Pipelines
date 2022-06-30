rule merge_units_R1:
    input:
        get_unit_R1_of_sample
    output:
        temp("results/merged_units/{sample}_R1.fastq.gz")
    params:
        files = lambda wc, input: " ".join(input)
    log:
        "logs/merged_units/{sample}_R1.log"
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell("cat {params.files} > results/merged_units/{wildcards.sample}_R1.fastq")
            shell("gzip results/merged_units/{wildcards.sample}_R1.fastq")

rule merge_units_R2:
    input:
        get_unit_R2_of_sample
    output:
        temp("results/merged_units/{sample}_R2.fastq.gz")
    params:
        files = lambda wc, input: " ".join(input)
    log:
        "logs/merged_units/{sample}_R2.log"
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell("cat {params.files} > results/merged_units/{wildcards.sample}_R2.fastq")
            shell("gzip results/merged_units/{wildcards.sample}_R2.fastq")

rule bwa_meth:
    input:
        reads = ["results/merged_units/{sample}_R1.fastq.gz"] if config["single_end"] else ["results/merged_units/{sample}_R1.fastq.gz","results/merged_units/{sample}_R2.fastq.gz"],
        idx = rules.bwa_index_meth.output
    output:
        temp("results/mapped/{sample}.sam")
    params:
        f"{assembly_path}{assembly}.fa"
    log:
        "logs/bwa/bwa_meth/{sample}.log"
    threads: 24
    conda:
        "../envs/bwa_meth.yaml"
    shell:
        "bwameth.py --threads {threads} --reference {params} {input.reads} > {output} 2>{log}"

rule sam_to_bam:
    input:
        "results/mapped/{sample}.sam",
    output:
        bam="results/mapped/{sample}.bam",
    log:
        "logs/bwa/{sample}_sam_to_bam.log",
    params:
        extra="-bS",  # optional params string
    threads: 2
    wrapper:
        "v1.7.0/bio/samtools/view"

rule samtools_sort_mapped:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.sorted.bam"
    params:
        extra=""
    log:
        "logs/bamtools_filtered/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule mark_merged_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/picard_dedup/{sample}.bam"),
        metrics="results/picard_dedup/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true PROGRAM_RECORD_ID='null' VALIDATION_STRINGENCY=LENIENT",
    wrapper:
        "v0.87.0/bio/picard/markduplicates"

rule bismark_map_pe:
    input:
        fq_1="results/merged_units/{sample}_R1.fastq.gz",
        fq_2="results/merged_units/{sample}_R2.fastq.gz",
        genome=f"{assembly_path}{assembly}.fa",
        bismark_indexes_dir=f"{assembly_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}_pe.bam"),
        report="results/bismark/reports/{sample}_PE_report.txt",
    log:
        "logs/bismark/map/{sample}.log"
    params:
        extra=f" --multicore 4 {config['params']['bismark']['map']['extra']}",
        basename='{sample}'
    threads: 24
    wrapper:
        "v1.7.0/bio/bismark/bismark"

rule bismark_map_se:
    input:
        fq="results/merged_units/{sample}_R1.fastq.gz",
        genome=f"{assembly_path}{assembly}.fa",
        bismark_indexes_dir=f"{assembly_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}_se.bam"),
        report="results/bismark/reports/{sample}_SE_report.txt",
    log:
        "logs/bismark/map/{sample}.log",
    params:
        # optional params string
        extra=f" --multicore 4 {config['params']['bismark']['map']['extra']}",
        basename='{sample}'
    threads: 24
    wrapper:
        "v1.7.0/bio/bismark/bismark"

rule deduplicate_bismark:
    input: 
        get_bismark_bam
    output:
        bam="results/bismark_mapped/{sample}.deduplicated.bam",
        report="results/bismark/reports/{sample}.deduplication_report.txt",
    log:
        "logs/bismark/{sample}.deduplicated.log",
    params:
        extra= lambda wc: "-s" if config["single_end"] else "-p"  # optional params string
    wrapper:
        "v1.7.0/bio/bismark/deduplicate_bismark"