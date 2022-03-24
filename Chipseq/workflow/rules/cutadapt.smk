rule cutadapt_pe:
    input:
        get_fastqs
    output:
        fastq1="results/trimmed/{sample}-{unit}_1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}_2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.pe.qc.txt",
        log="logs/cutadapt/{sample}-{unit}.pe.log"
    params:
        adapters = config["params"]["cutadapt-pe"],
        extra = config["params"]["cutadapt-others"]
    log:
        "logs/cutadapt/{sample}-{unit}.pe.log"
    wrapper:
        "v1.3.1/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_fastqs
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.se.qc.txt",
        log="logs/cutadapt/{sample}-{unit}.se.log"
    params:
        adapters="{}".format(
            config["params"]["cutadapt-se"]),
        extra = "{}".format(
            config["params"]["cutadapt-others"]
            )
    log:
        "logs/cutadapt/{sample}-{unit}.se.log"
    wrapper:
        "v1.3.1/bio/cutadapt/se"