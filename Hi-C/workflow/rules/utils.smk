rule samtools_index:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step}/{samples_units}.bam.bai"
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/{step}/{samples_units}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"

rule fanc_fragments:
    input:
        f"{assembly_path}{assembly}.fa"
    output:
        f"resources/{assembly}.{enzyme_file}.{fragments_file}.fragments.bed"
    params:
        chr = f"-c {config['params']['fanc']['chr']}" if config['params']['fanc']['chr'] else "",
        enzyme = config["params"]["fanc"]["enzyme"]
    log:
        f"logs/ref/{assembly}.{enzyme_file}.fragments.log"
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc fragments {params.chr} {input} {params.enzyme} {output} 2>{log}"