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

rule bowtie2:
    input:
        sample=["results/merged_units/{sample}_R1.fastq.gz"] if config["single_end"] else ["results/merged_units/{sample}_R1.fastq.gz","results/merged_units/{sample}_R2.fastq.gz"],
        idx=multiext(
            f"{assembly_path}{assembly}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("results/mapped/{sample}.bam")
    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra="--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700",  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "v1.5.0/bio/bowtie2/align"

rule bowtie2_spike:
    input:
        sample=["results/merged_units/{sample}_R1.fastq.gz"] if config["single_end"] else ["results/merged_units/{sample}_R1.fastq.gz","results/merged_units/{sample}_R2.fastq.gz"],
        idx=multiext(
            f"{spike_path}{spike_assembly}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("results/mapped/{sample}_spike-in.bam")
    log:
        "logs/bowtie2/{sample}_spike-in.log",
    params:
        extra="--no-overlap --no-dovetail --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700",  # optional parameters
    threads: 8  # Use at least two threads
    conda:
        "../envs/bowtie2.yaml"
    script:
        "../scripts/bowtie2.py"

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
        bams="results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/picard_dedup/{sample}.sorted.bam"),
        metrics="results/picard_dedup/{sample}.sorted.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        extra=f"""--REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT {f"--OPTICAL_DUPLICATE_PIXEL_DISTANCE {config['params']['optical_distance']}" if config['params']['optical_distance'] > 0 else ''}""",
    resources:
        mem_mb=16384,
    threads: 4
    wrapper:
        "v3.0.2/bio/picard/markduplicates"
