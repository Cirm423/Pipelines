rule sort_mapped:
    input:
        "results/mapped/{sample}_R{read}.bam"
    output:
        "results/mapped/{sample}_R{read}.sorted.bam"
    params:
        extra="-n"
    log:
        "logs/bamtools_filtered/{sample}_R{read}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule fanc_pairs:
    input:
        R1 = "results/mapped/{sample}_R1.sorted.bam",
        R2 = "results/mapped/{sample}_R2.sorted.bam",
        genome = f"resources/{assembly}.{enzyme_file}.fragments.bed"
    output:
        pairs = "results/pairs/{sample}.pairs",
        stats = report("results/pairs/{sample}.pairs_stats.pdf",category="Pairs"),
        dist_plot = report("results/pairs/{sample}.re-dist.png",category="Pairs"),
        l_error = report("results/pairs/{sample}.ligation-err.png",category="Pairs")
    params:
        unmap = "-m" if config["params"]["fanc"]["filter"]["unmap"] else "",
        multimap = "" if not config["params"]["fanc"]["filter"]["multimap"] else ("-us" if config["params"]["fanc"]["filter"]["multimap"] == "strict" else "-u"),
        inward = f"-i {config['params']['fanc']['filter']['inward']}" if config['params']['fanc']['filter']['inward'] else "",
        outward = f"-o {config['params']['fanc']['filter']['outward']}" if config['params']['fanc']['filter']['outward'] else "",
        distance = f"-d {config['params']['fanc']['filter']['distance']}" if config['params']['fanc']['filter']['distance'] else "",
        ligation = f"-l {config['params']['fanc']['filter']['ligation']}" if config['params']['fanc']['filter']['ligation'] else "",
        pcr = f"-p {config['params']['fanc']['filter']['pcr']}" if config['params']['fanc']['filter']['pcr'] else "",
        extra = config['params']['fanc']['filter']['extra']
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc pairs {input.R1} {input.R2} {input.genome} \
        {params.unmap} {params.multimap} {params.inward} {params.outward} \
        {params.distance} {params.ligation} {params.pcr} {params.extra} \
        --statistics-plot {output.stats} --re-dist-plot {output.dist_plot} --ligation-error-plot {output.l_error}"""

rule fanc_hic:
    input:
        "results/pairs/{sample}.pairs", #Maybe change to groups to account for biological replicates
    output:
        hic = "results/hic/{sample}.hic",
        stats = report("results/hic/{sample}.hic_stats.pdf",category="Hi-C matrix")
    params:
        bin = f"-b {config['params']['fanc']['hic']['bin_size']}",
        filter = config['params']['fanc']['hic']['filter'],
        diag = f"-d {config['params']['fanc']['hic']['diagonal']}" if config['params']['fanc']['hic']['bin_size'] else "",
        norm = (f"-n -m {config['params']['fanc']['hic']['normalize']['method']} -w" if config['params']['fanc']['hic']['normalize']['whole'] else f"-n -m {config['params']['fanc']['hic']['normalize']['method']}") if config['params']['fanc']['hic']['normalize']['activate'] else "",
        extra = config['params']['fanc']['hic']['extra']
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc hic {input} {output.hic} --statistics-plot {output.stats} \
        {params.bin} {params.filter} {params.diag} {params.norm} {params.extra} -t {threads}"""