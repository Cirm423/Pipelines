rule methyldackel:
    input:
        bam = "results/picard_dedup/{sample}.bam",
        bai = "results/picard_dedup/{sample}.bam.bai",
        fa = f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        fai = f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa.fai",
    output:
        bg = "results/methyldackel/{sample}_CpG.bedgraph",
        txt = "results/methyldackel/{sample}_methyldackel.txt"
        svg = report("results/methyldackel/{sample}_meth_bias.svg", category="Methylation")
    params:
        prefix = lambda wc: "results/methyldackel/{wc.sample}",
        prefix_mbias = lambda wc: "results/methyldackel/{wc.sample}_meth_bias",
        comprehensive = "--CHG --CHH" if config["params"]["methyldackel"]["comprehensive"] else "",
        min_depth = config["params"]["methyldackel"]["min_depth"] if config["params"]["methyldackel"]["min_depth"] > 0 else "",
        ignore_flags = "--ignoreFlags" if config["params"]["methyldackel"]["ignore_flags"] else "",
        methyl_kit = "--methylKit" if config["params"]["methyldackel"]["methyl_kit"] else "",
        extra_extract = config["params"]["methyldackel"]["extra"]["extract"],
        extra_mbias = config["params"]["methyldackel"]["extra"]["mbias"],
    conda:
        "../envs/methyldackel.yaml"
    threads: 24
    shell:
        "MethyDackel extract -@ {threads} {params.comprehensive} {params.ignore_flags} {params.methyl_kit} {params.min_depth} {input.fa} {input.bam} {params.extra_extract} -o {params.prefix}"
        "MethyDackel mbias -@ {threads} {params.comprehensive} {params.ignore_flags} {params.extra_mbias} {input.fa} {input.bam} {params.prefix_mbias} --txt > {output.txt}"

