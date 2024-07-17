## RSEQC

rule rseqc_gtf2bed:
    input:
        f"{assembly_path}{assembly}.annotation.gtf",
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"

rule rseqc_junction_annotation:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bed="results/qc/rseqc/annotation.bed",
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end)
    output:
        "results/qc/rseqc/{sample}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bed="results/qc/rseqc/annotation.bed",
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    threads: 24
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}.log",
    conda:
        "../envs/rseqc2.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bed="results/qc/rseqc/annotation.bed",
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}.log",
    conda:
        "../envs/rseqc2.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bed="results/qc/rseqc/annotation.bed",
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc3.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bed="results/qc/rseqc/annotation.bed",
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}.log",
    conda:
        "../envs/rseqc3.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


# rule rseqc_readdup:
#     input:
#         bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
#         bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
#     output:
#         "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
#     priority: 1
#     log:
#         "logs/rseqc/rseqc_readdup/{sample}.log",
#     params:
#         prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
#     conda:
#         "../envs/rseqc.yaml"
#     shell:
#         "read_duplication.py -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        bam=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam",lib_end=lib_end),
        bai=expand("results/star/{lib_end}/{{sample}}/Aligned.sortedByCoord.out.bam.bai",lib_end=lib_end),
    output:
        "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc4.yaml"
    threads: 24
    shell:
        "read_GC.py -i {input.bam} -o {params.prefix} > {log} 2>&1"

rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    params:
        extra = "--quiet"
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
        extra = "--quiet"
    log:
        "logs/fastqc/trimmed_{sample}.{unit}.{read}.log"
    threads: 6
    resources:
        mem_mb = 2048
    wrapper:
        "v2.6.0/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input,
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/multiqc.py"
