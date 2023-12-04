rule get_sra_pe:
    output:
        "sra-pe-reads/{accession}_1.fastq.gz",
        "sra-pe-reads/{accession}_2.fastq.gz",
    params:
        extra="--skip-technical",
    log:
        "logs/get-sra/{accession}.log",
    threads: 6 
    wrapper:
        "v2.6.0/bio/sra-tools/fasterq-dump"


rule get_sra_se:
    output:
        "sra-se-reads/{accession}.fastq.gz",
    params:
        extra="--skip-technical",
    log:
        "logs/get-sra/{accession}.log",
    threads: 6 
    wrapper:
        "v2.6.0/bio/sra-tools/fasterq-dump"

rule fix_sra_se:
    input:
        "sra-se-reads/{accession}.fastq.gz",
    output:
        "sra-se-reads/{accession}.fixed.fastq.gz",
    log:
        "logs/fix-sra-se/{accession}.log"
    threads: 6
    shell:
        "zcat {input} | sed -e 's/^\(@[^[:blank:]]*\)[[:blank:]]\+/\1_/' -e 's/^\(+[^[:blank:]]*\)[[:blank:]]\+/\1_/' | gzip > {output}"

rule fix_sra_pe:
    input:
        read1="sra-pe-reads/{accession}_1.fastq.gz",
        read2="sra-pe-reads/{accession}_2.fastq.gz"
    output:
        read1="sra-pe-reads/{accession}_1.fixed.fastq.gz",
        read2="sra-pe-reads/{accession}_2.fixed.fastq.gz"
    log:
        "logs/fix-sra-pe/{accession}.log"
    threads: 6
    shell:
        """zcat {input.read1} | sed -e 's/^\(@[^[:blank:]]*\)[[:blank:]]\+/\1_/' -e 's/^\(+[^[:blank:]]*\)[[:blank:]]\+/\1_/' | gzip > {output.read1} && \
        zcat {input.read2} | sed -e 's/^\(@[^[:blank:]]*\)[[:blank:]]\+/\1_/' -e 's/^\(+[^[:blank:]]*\)[[:blank:]]\+/\1_/' | gzip > {output.read2}"""

rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        temp("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
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
        extra=config["params"]["cutadapt-others"],
        adapters=config["params"]["cutadapt-pe"],
    threads: 8
    wrapper:
        "v2.6.0/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed/{sample}_{unit}_single.fastq.gz",
        qc="results/trimmed/{sample}_{unit}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra=config["params"]["cutadapt-others"],
        adapters_r1=config["params"]["cutadapt-se"],
    threads: 8
    wrapper:
        "v2.6.0/bio/cutadapt/se"


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