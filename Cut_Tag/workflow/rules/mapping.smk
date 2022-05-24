rule merge_units_R1:
    input:
        get_unit_R1_of_sample
    output:
        "results/merged_units/{sample}_R1.fastq.gz"
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
        "results/merged_units/{sample}_R2.fastq.gz"
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
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}",
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
            f"{config['resources']['path']}{config['resources']['ref']['spike_assembly']}",
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
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "v1.5.0/bio/bowtie2/align"

rule mark_merged_duplicates:
    input:
        "results/merged/{sample}.bam"
    output:
        bam=temp("results/picard_dedup/{sample}.bam"),
        metrics="results/picard_dedup/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT",
    wrapper:
        "v0.87.0/bio/picard/markduplicates"
