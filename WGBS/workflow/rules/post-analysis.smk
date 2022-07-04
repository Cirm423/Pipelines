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
    threads: 8
    conda:
        "../envs/preseq.yaml"
    script:
        "../scripts/preseq.py"