rule format_to_sam:
    input:
        get_deduplicated_bams
    output:
        temp("results/methpipe/sams/{sample}.deduplicated.sam")
    params:
        "-f bismark" if config["params"]["mode"] == "bismark" else "-f general"
    log:
        "logs/methpipe/{sample}.format_to_sam.log"
    threads: 24
    conda:
        "../envs/methpipe.yaml"
    shell:
        "format_reads {params} -o {output} {input} 2>{log}" 

rule samtools_sort_sams:
    input:
        "results/methpipe/sams/{sample}.deduplicated.sam"
    output:
        temp("results/methpipe/sams/{sample}.deduplicated.sorted.sam")
    params:
        extra=""
    log:
        "logs/methpipe/{sample}_sort_sam.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule meth_counts:
    input:
        sam="results/methpipe/sams/{sample}.deduplicated.sorted.sam",
        ref=f"{assembly_path}{assembly}.fa"
    output:
        temp("results/methpipe/meth_counts/{sample}.meth")
    log:
        "logs/methpipe/{sample}.meth_counts.log"
    threads: 24
    conda:
        "../envs/methpipe.yaml"
    shell:
        "methcounts -c {input.ref} -n -o {output} {input.sam} 2>{log}"

# rule sort_meth_counts:
#     input:
#         "results/methpipe/meth_counts/{sample}.meth"
#     output:
#         "results/methpipe/meth_counts/{sample}.sorted.meth"
#     log:
#         "logs/methpipe/{sample}.sort_counts.log"
#     threads: 8
#     conda:
#         "../envs/coreutils.yaml"
#     shell:
#         "sort --parallel={threads} -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o {output} {input} 2>{log}"

rule symmetric_cpgs:
    input:
        "results/methpipe/meth_counts/{sample}.meth"
    output:
        "results/methpipe/meth_counts/{sample}.symmetric.meth"
    log:
        "logs/methpipe/{sample}.symmetric.log"
    threads: 24
    conda:
        "../envs/methpipe.yaml"
    shell:
        "symmetric-cpgs -o {output} {input}"

rule call_PMDs:
    input:
        "results/methpipe/meth_counts/{sample}.symmetric.meth"
    output:
        "results/methpipe/PMDs/{sample}.pmd"
    log:
        "logs/methpipe/{sample}.PMDs.log"
    threads: 24
    conda:
        "../envs/methpipe.yaml"
    shell:
        "pmd -i 1000 -o {output} {input} 2>{log}"