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
    wrapper:
        "v1.3.1/bio/fastqc"

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
    wrapper:
        "v1.3.1/bio/fastqc"

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

rule ATACseqQC:
    input:
        "results/filtered/{sample}.sorted.bam"
    output:
        fragmentSizeDistribution = report("results/qc/ATACseqQC/fragmentSizeDistribution.pdf", category = "ATACseqQC"),
        
    params:
        BSgenome = "",
        Txdb = ""
    log:
        "logs/ATACseqQC/{sample}.log"
    conda:
        "../envs/ATACseqQC.yaml"
    script:
        "../scripts/ATACseqQC.R"