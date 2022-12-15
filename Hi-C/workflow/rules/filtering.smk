rule sort_mapped:
    input:
        "results/mapped/{sample}_R{read}.bam"
    output:
        "results/mapped/{sample}_R{read}.sorted.bam"
    params:
        extra="-n"
    log:
        "logs/bamtools_filtered/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"
        
rule fanc_pairs:
    input:
        R1 = "results/mapped/{sample}_R1.sorted.bam",
        R2 = "results/mapped/{sample}_R2.sorted.bam",
        genome = f"{assembly_path}{assembly}.fa"
    output:
        pairs = "results/pairs/{sample}.pairs",
        stats = report("results/pairs/{sample}.stats.pdf",category="Filter")
    params:
        unmap = "-m" if config["params"]["fanc"]["filter"]["unmap"] else "",
        multimap = "" if not config["params"]["fanc"]["filter"]["multimap"] else ("-us" if config["params"]["fanc"]["filter"]["multimap"] == "strict" else "-u"),
        inward = f"-i {config['params']['fanc']['filter']['inward']}" if config['params']['fanc']['filter']['inward'] else "",
        outward = f"-o {config['params']['fanc']['filter']['outward']}" if config['params']['fanc']['filter']['outward'] else "",
        distance = f"-d {config['params']['fanc']['filter']['distance']}" if config['params']['fanc']['filter']['distance'] else "",
        ligation = f"-l {config['params']['fanc']['filter']['ligation']}" if config['params']['fanc']['filter']['ligation'] else "",
        pcr = f"-p {config['params']['fanc']['filter']['pcr']}" if config['params']['fanc']['filter']['pcr'] else ""
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc pairs {input.R1} {input.R2} {input.genome} \
        {params.unmap} {params.multimap} {params.inward} {params.outward} \
        {params.distance} {params.ligation} {params.pcr} \
        --statistics-plot {output.stats} --re-dist-plot --ligation-error-plot"""