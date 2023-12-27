rule merge_units_R1:
    input:
        get_unit_R1_of_sample,
    output:
        temp("results/merged_units/{sample}_R1.fastq.gz"),
    params:
        files=lambda wc, input: " ".join(input),
    log:
        "logs/merged_units/{sample}_R1.log",
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell(
                "cat {params.files} > results/merged_units/{wildcards.sample}_R1.fastq"
            )
            shell("gzip results/merged_units/{wildcards.sample}_R1.fastq")


rule merge_units_R2:
    input:
        get_unit_R2_of_sample,
    output:
        temp("results/merged_units/{sample}_R2.fastq.gz"),
    params:
        files=lambda wc, input: " ".join(input),
    log:
        "logs/merged_units/{sample}_R2.log",
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell(
                "cat {params.files} > results/merged_units/{wildcards.sample}_R2.fastq"
            )
            shell("gzip results/merged_units/{wildcards.sample}_R2.fastq")


rule bwa_meth:
    input:
        reads=["results/merged_units/{sample}_R1.fastq.gz"] if config["single_end"] else ["results/merged_units/{sample}_R1.fastq.gz","results/merged_units/{sample}_R2.fastq.gz"],
        idx=rules.bwa_index_meth.output,
    output:
        temp("results/mapped/{sample}.sam"),
    params:
        f"{assembly_path}{assembly}.fa",
    log:
        "logs/bwa/bwa_meth/{sample}.log",
    threads: 24
    conda:
        "../envs/bwa_meth.yaml"
    shell:
        "bwameth.py --threads {threads} --reference {params} {input.reads} > {output} 2>{log}"


rule samtools_sort_mapped:
    input:
        "results/mapped/{sample}.sam",
    output:
        temp("results/mapped/{sample}.bam"),
    params:
        extra="",
    log:
        "logs/mapped/{sample}_sam_to_sbam.log",
    threads: 8
    wrapper:
        "v1.3.1/bio/samtools/sort"


rule mark_merged_duplicates:
    input:
        bams="results/mapped/{sample}.bam",
    output:
        bam="results/picard_dedup/{sample}.bam",
        metrics="results/picard_dedup/{sample}.metrics.txt",
    log:
        "logs/picard/picard_dedup/{sample}.log",
    params:
        extra=f"--REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE {config['params']['methyldackel']['optical_distance']}",
    threads: 24
    wrapper:
        "v3.0.2/bio/picard/markduplicates"


rule bismark_map_pe:
    input:
        fq_1="results/merged_units/{sample}_R1.fastq.gz",
        fq_2="results/merged_units/{sample}_R2.fastq.gz",
        genome=f"{assembly_path}{assembly}.fa",
        bismark_indexes_dir=f"{assembly_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}_pe.bam"),
        report="results/bismark_mapped/{sample}_PE_report.txt",
        tmp_dir=temp(directory("tmp/{sample}_bmap/"))
    log:
        "logs/bismark/map/{sample}_pe.log",
    params:
        extra=f" {config['params']['bismark']['map']['extra']} --temp_dir tmp/{{sample}}_bmap/",
        #basename='{sample}'
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/bismark"


rule bismark_map_se:
    input:
        fq="results/merged_units/{sample}_R1.fastq.gz",
        genome=f"{assembly_path}{assembly}.fa",
        bismark_indexes_dir=f"{assembly_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}_se.bam"),
        report="results/bismark_mapped/{sample}_SE_report.txt",
        tmp_dir=temp(directory("tmp/{sample}_bmap/"))
    log:
        "logs/bismark/map/{sample}_se.log",
    params:
        # optional params string
        extra=f" {config['params']['bismark']['map']['extra']} --temp_dir tmp/{{sample}}_bmap/",
        #basename='{sample}'
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/bismark"


rule deduplicate_bismark_pe:
    input:
        "results/bismark_mapped/{sample}_pe.bam",
    output:
        bam="results/bismark_mapped/{sample}_pe.deduplicated.bam",
        report="results/bismark_mapped/{sample}_pe.deduplication_report.txt",
    log:
        "logs/bismark/{sample}_pe.deduplicated.log",
    params:
        extra="-p",  # optional params string
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/deduplicate_bismark"


rule deduplicate_bismark_se:
    input:
        "results/bismark_mapped/{sample}_se.bam",
    output:
        bam="results/bismark_mapped/{sample}_se.deduplicated.bam",
        report="results/bismark_mapped/{sample}.deduplication_report.txt",
    log:
        "logs/bismark/{sample}_se.deduplicated.log",
    params:
        extra="-s",  # optional params string
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/deduplicate_bismark"


rule bismark_map_phage_pe:
    input:
        fq_1="results/merged_units/{sample}_R1.fastq.gz",
        fq_2="results/merged_units/{sample}_R2.fastq.gz",
        genome=f"{phage_path}phage_lambda.fa",
        bismark_indexes_dir=f"{phage_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}-phage_pe.bam"),
        report="results/bismark_mapped/{sample}-phage_PE_report.txt",
        tmp_dir=temp(directory("tmp/{sample}_bmap_phage/"))
    log:
        "logs/bismark/map/{sample}-phage_pe.log",
    params:
        extra=f" {config['params']['bismark']['map']['extra']} --temp_dir tmp/{{sample}}_bmap_phage/",
        #basename="{sample}",
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/bismark"


rule bismark_map_phage_se:
    input:
        fq="results/merged_units/{sample}_R1.fastq.gz",
        genome=f"{phage_path}phage_lambda.fa",
        bismark_indexes_dir=f"{phage_path}Bisulfite_Genome",
        #genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam=temp("results/bismark_mapped/{sample}-phage_se.bam"),
        report="results/bismark_mapped/{sample}-phage_SE_report.txt",
        tmp_dir=temp(directory("tmp/{sample}_bmap_phage/"))
    log:
        "logs/bismark/map/{sample}-phage_se.log",
    params:
        # optional params string
        extra=f" {config['params']['bismark']['map']['extra']} --temp_dir tmp/{{sample}}_bmap_phage/",
        basename="{sample}",
    threads: 24
    wrapper:
        "v2.2.0/bio/bismark/bismark"
