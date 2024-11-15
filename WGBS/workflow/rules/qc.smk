rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    threads: 6
    resources:
        mem_mb = 2048
    wrapper:
        "v2.6.0/bio/fastqc"

rule fastqc_trimmed:
    input:
        get_individual_trimmed_fastq
    output:
        html="results/qc/fastqc/trimmed_{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/trimmed_{sample}.{unit}.{read}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/trimmed_{sample}.{unit}.{read}.log"
    threads: 6
    resources:
        mem_mb = 2048
    wrapper:
        "v2.6.0/bio/fastqc"

rule sort_deduplicated:
    input:
        get_dedup_bam
    output:
        temp("results/qualimap/{sample}.sorted.bam")
    params:
        extra=""
    log:
        "logs/qualimap/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule qualimap:
    input:
        "results/qualimap/{sample}.sorted.bam"
    output:
        directory("results/qualimap/{sample}_qualimap")
    log:
        "logs/qualimap/{sample}.log"
    params:
        gcref = "-gd HUMAN" if "h" in assembly else "-gd MOUSE" if "m" in assembly else ""
    threads: 8
    conda:
        "../envs/qualimap.yaml"
    shell:
        "JAVA_OPTS='-Djava.awt.headless=true' qualimap bamqc {params.gcref} -bam {input} -outdir {output} --collect-overlap-pairs -nt {threads} --java-mem-size=16G"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "results/qc/multiqc/multiqc.html"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/multiqc.py"