## RSEQC

rule rseqc_gtf2bed:
    input:
        f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
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
        bam=lambda wc: get_star_bam(wc,original=True),
        bed="results/qc/rseqc/annotation.bed",
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.junctionanno.junction.bed"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_junction_annotation/?.log"),
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
        bam=lambda wc: get_star_bam(wc,original=True),
        bed="results/qc/rseqc/annotation.bed",
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.junctionsat.junctionSaturation_plot.pdf"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_junction_saturation/?.log"),
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
        bam=lambda wc: get_star_bam(wc,original=True),
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.stats.txt"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_stat/?.log"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam=lambda wc: get_star_bam(wc,original=True),
        bed="results/qc/rseqc/annotation.bed",
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.infer_experiment.txt"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_infer/?.log"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=lambda wc: get_star_bam(wc,original=True),
        bed="results/qc/rseqc/annotation.bed",
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.inner_distance_freq.inner_distance.txt"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_innerdis/?.log"),
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=lambda wc: get_star_bam(wc,original=True),
        bed="results/qc/rseqc/annotation.bed",
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.readdistribution.txt"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_readdis/?.log"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


# rule rseqc_readdup:
#     input:
#         bam=lambda wc: get_star_bam(wc,original=True),
#         bai=lambda wc: get_star_bam_bai(wc,original=True),
#     output:
#         path_merged_cond("results/qc/rseqc/?.readdup.DupRate_plot.pdf"),
#     priority: 1
#     log:
#         path_merged_cond("logs/rseqc/rseqc_readdup/?.log"),
#     params:
#         prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
#     conda:
#         "../envs/rseqc.yaml"
#     shell:
#         "read_duplication.py -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        bam=lambda wc: get_star_bam(wc,original=True),
        bai=lambda wc: get_star_bam_bai(wc,original=True),
    output:
        path_merged_cond("results/qc/rseqc/?.readgc.GC_plot.pdf"),
    priority: 1
    log:
        path_merged_cond("logs/rseqc/rseqc_readgc/?.log"),
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
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
        ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    threads: 6
    wrapper:
        "0.72.0/bio/fastqc"

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
        "0.72.0/bio/fastqc"

rule multiqc:
    input:
        lambda wc: get_star_output_all_units(wc, fi="bam", orig =True),
        get_multiqc_input,
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.junctionanno.junction.bed"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.junctionsat.junctionSaturation_plot.pdf"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.infer_experiment.txt"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.stats.txt"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.inner_distance_freq.inner_distance.txt"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.readdistribution.txt"),
            unit=units.itertuples(),
        ))),
        # list(set(expand(
        #     path_merged_cond_mqc("results/qc/rseqc/?.readdup.DupRate_plot.pdf"),
        #     unit=units.itertuples(),
        # ))),
        list(set(expand(
            path_merged_cond_mqc("results/qc/rseqc/?.readgc.GC_plot.pdf"),
            unit=units.itertuples(),
        ))),
        list(set(expand(
            path_merged_cond_mqc("logs/rseqc/rseqc_junction_annotation/?.log"),
            unit=units.itertuples(),
        ))),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/multiqc.py"
