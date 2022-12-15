rule merge_units_R1:
    input:
        get_unit_R1_of_sample
    output:
        temp("results/merged_units/{sample}_R1.fastq.gz")
    params:
        files = lambda wc, input: " ".join(input)
    log:
        "logs/merged_units/{sample}_R1.log"
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell("cat {params.files} > results/merged_units/{wildcards.sample}_R1.fastq")
            shell("gzip results/merged_units/{wildcards.sample}_R1.fastq")

rule merge_units_R2:
    input:
        get_unit_R2_of_sample
    output:
        temp("results/merged_units/{sample}_R2.fastq.gz")
    params:
        files = lambda wc, input: " ".join(input)
    log:
        "logs/merged_units/{sample}_R2.log"
    run:
        if input[0].endswith("gz"):
            shell("cat {params.files} > {output} 2>{log}")
        else:
            shell("cat {params.files} > results/merged_units/{wildcards.sample}_R2.fastq")
            shell("gzip results/merged_units/{wildcards.sample}_R2.fastq")

rule fanc_map:
    input:
        read = "results/merged_units/{sample}_R{read}.fastq.gz",
        idx = rules.bwa_index.output
    output:
        temp("results/mapped/{sample}_R{read}.bam")
    params:
        enzyme = onfig["params"]["fanc"]["enzyme"],
        min_size = config["params"]["fanc"]["map"]["min_size"],
        step_size = config["params"]["fanc"]["map"]["step_size"],
        quality = config["params"]["fanc"]["map"]["quality"],
        extra = config["params"]["fanc"]["map"]["extra"],
    log:
        "logs/map/{sample}_R{read}.log"
    threads: 16
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc map {input.read} {input.idx} {output} \
        -r {params.enzyme} -m {params.min_size} -s {params.step_size} -q {params.quality} \
        {params.extra} -t {threads} --mapper-type bwa"""