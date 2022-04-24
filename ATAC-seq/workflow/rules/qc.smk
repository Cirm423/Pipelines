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
        LibComplexity = report("results/qc/ATACseqQC/{sample}/LibComplexity.pdf", category = "ATACseqQC_{sample}"),
        fragmentSizeDistribution = report("results/qc/ATACseqQC/{sample}/fragmentSizeDistribution.pdf", category = "ATACseqQC_{sample}"),
        PTscore = report("results/qc/ATACseqQC/{sample}/PTscore.pdf", category = "ATACseqQC_{sample}"),
        NFRscore = report("results/qc/ATACseqQC/{sample}/NFRscore.pdf", category = "ATACseqQC_{sample}"),
        TSSEscore = report("results/qc/ATACseqQC/{sample}/TSSEscore.pdf", category = "ATACseqQC_{sample}"),
        cumulativePercentage = report("results/qc/ATACseqQC/{sample}/cumulativePercentage.pdf", category = "ATACseqQC_{sample}"),
        featureAligndHeatmap = report("results/qc/ATACseqQC/{sample}/featureAligndHeatmap.pdf", category = "ATACseqQC_{sample}"),
        TSS_profile = report("results/qc/ATACseqQC/{sample}/TSS_profile.pdf", category = "ATACseqQC_{sample}"),
        CTCF_footprint = report("results/qc/ATACseqQC/{sample}/CTCF_footprint.pdf", category = "ATACseqQC_{sample}"),
        CTCF_Vplot = report("results/qc/ATACseqQC/{sample}/CTCF_Vplot.pdf", category = "ATACseqQC_{sample}")
    params:
        BSgenome = QC_packages[config['resources']['ref']['assembly']]["BSgenome"],
        Txdb = QC_packages[config['resources']['ref']['assembly']]["Txdb"]
    log:
        "logs/ATACseqQC/{sample}.log"
    conda:
        "../envs/ATACseqQC.yaml"
    script:
        "../scripts/ATACseqQC.R"