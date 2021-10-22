rule get_sra_pe:
    output:
        "sra-pe-reads/{accession}_1.fastq.gz",
        "sra-pe-reads/{accession}_2.fastq.gz",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"


rule get_sra_se:
    output:
        "sra-se-reads/{accession}.fastq.gz",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"

# rule compress_sra:
#     input:
#         "sra/{accession}_{read}.fastq",
#     output:
#         "sra/{accession}_{read}.fastq.gz",
#     log:
#         "logs/compress-sra/{accession}_{read}.log"
#     shell:
#         "gzip {input} 2> {log}"

rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}_{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}_{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        others=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed/{sample}_{unit}_single.fastq.gz",
        qc="results/trimmed/{sample}_{unit}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        others=config["params"]["cutadapt-se"],
        adapters_r1=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/se"


rule merge_fastqs_gz:
    input:
        get_fastqs_gz,
    output:
        "results/merged/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"

rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "results/merged/{sample}_{read}.fastq",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"