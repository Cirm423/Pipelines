rule samtools_index:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        temp("results/{step}/{samples_units}.bam.bai")
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/{step}/{samples_units}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"

rule fanc_fragments:
    input:
        fasta = f"{assembly_path}{assembly}.fa",
        ini = "results/fanc/package.done"
    output:
        "resources/{assembly}.{enzyme}.{fragments}.fragments.bed"
    params:
        chr = f"-c {fragment_chr}" if fragment_chr else "",
        enzyme = config["params"]["fanc"]["enzyme"]
    log:
        "logs/ref/{assembly}.{enzyme}.{fragments}.log"
    threads: 4
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc fragments {params.chr} {input.fasta} {params.enzyme} {output} 2>{log}"