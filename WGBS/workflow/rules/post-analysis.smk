
rule samtools_sort_picard_dedup:
    input:
        "results/picard_dedup/{sample}.bam"
    output:
        "results/picard_dedup/{sample}.sorted.bam"
    params:
        extra=""
    log:
        "logs/bamtools_filtered/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule qualimap:
    input:
        get_dedup_bam
    output:
        "results/qc/qualimap/{sample}"
    params:
        genome = lambda wc: "-gd HUMAN" if "h" in assembly else "-gd MOUSE" if "m" in assembly else "",
    log:
        "logs/qc/{sample}_qualimap.log"
    conda:
        "../envs/qualimap.yaml"
    threads: 8
    shell:
        "qualimap bamqc {params.genome} -bam {input} -outdir {output} --collect-overlap-pairs -nt {threads}"

rule preseq_lc_extrap:
    input:
        get_dedup_bam
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        lambda wildcards, resources: f"-v {'' if config['single_end'] else '-pe -seg_len 1000000000'} -seed 1 {'-D' if resources.attempt > 1 else ''}" 
    log:
        "logs/preseq/{sample}.log"
    resources:
        attempt = lambda wildcards, attempt: attempt
    conda:
        "../envs/preseq.yaml"
    script:
        "../scripts/preseq.py"