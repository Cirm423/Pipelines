rule sort_mapped:
    input:
        "results/mapped/{sample}_R{read}.bam"
    output:
        "results/mapped/{sample}_R{read}.{enzyme}.sorted.bam"
    params:
        extra="-n"
    log:
        "logs/bamtools_filtered/{sample}_R{read}_{enzyme}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule fanc_pairs:
    input:
        R1 = "results/mapped/{sample}_R1.{enzyme}.sorted.bam",
        R2 = "results/mapped/{sample}_R2.{enzyme}.sorted.bam",
        genome = f"resources/{assembly}.{enzyme_file}.{fragments_file}.fragments.bed"
    output:
        pairs = "results/pairs/{sample}.{enzyme}.{fragments}.pairs",
        stats = report("results/pairs/{sample}.{enzyme}.{fragments}.pairs_stats.pdf",category="Pairs"),
        dist_plot = report("results/pairs/{sample}.{enzyme}.{fragments}.re-dist.png",category="Pairs"),
        l_error = report("results/pairs/{sample}.{enzyme}.{fragments}.ligation-err.png",category="Pairs")
    params:
        unmap = "-m" if config["params"]["fanc"]["filter"]["unmap"] else "",
        multimap = "" if not config["params"]["fanc"]["filter"]["multimap"] else ("-us" if config["params"]["fanc"]["filter"]["multimap"] == "strict" else "-u"),
        inward = f"-i {config['params']['fanc']['filter']['inward']}" if config['params']['fanc']['filter']['inward'] else "",
        outward = f"-o {config['params']['fanc']['filter']['outward']}" if config['params']['fanc']['filter']['outward'] else "",
        distance = f"-d {config['params']['fanc']['filter']['distance']}" if config['params']['fanc']['filter']['distance'] else "",
        ligation = f"-l" if config['params']['fanc']['filter']['ligation'] else "",
        pcr = f"-p {config['params']['fanc']['filter']['pcr']}" if config['params']['fanc']['filter']['pcr'] else "",
        extra = config['params']['fanc']['filter']['extra']
    log:
        "logs/fanc/{sample}.{enzyme}.{fragments}.pairs.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc pairs {input.R1} {input.R2} {output.pairs} -g {input.genome} \
        {params.unmap} {params.multimap} {params.inward} {params.outward} \
        {params.distance} {params.ligation} {params.pcr} {params.extra} \
        --statistics-plot {output.stats} --re-dist-plot {output.dist_plot} --ligation-error-plot {output.l_error} 2>{log}"""

rule fanc_hic:
    input:
        get_pairs_files,
    output:
        "results/hic/{sample_group}.{enzyme}.{fragments}.fragment_level.hic",
    log:
        "logs/fanc/{sample_group}.{enzyme}.{fragments}.fragment_level.hic.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc hic {input} {output} -t {threads} 2>{log}"

rule fanc_hic_bin:
    input:
        "results/hic/{sample_group}.{enzyme}.{fragments}.fragment_level.hic",
    output:
        hic = "results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic",
        stats = report("results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic_stats.pdf",category="Hi-C matrix")
    params:
        bin = lambda wildcards: f"-b {wildcards.resolution}",
        filter = config['params']['fanc']['hic']['filter'],
        diag = f"-d {config['params']['fanc']['hic']['diagonal']}" if config['params']['fanc']['hic']['diagonal'] else "",
        norm = (f"-n -m {config['params']['fanc']['hic']['normalize']['method']} -w" if config['params']['fanc']['hic']['normalize']['whole'] else f"-n -m {config['params']['fanc']['hic']['normalize']['method']}") if config['params']['fanc']['hic']['normalize']['activate'] else "",
        extra = config['params']['fanc']['hic']['extra']
    log:
        "logs/fanc/{sample_group}.{enzyme}.{fragments}-{resolution}.hic.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        """fanc hic {input} {output.hic} --statistics-plot {output.stats} \
        {params.bin} {params.filter} {params.diag} {params.norm} {params.extra} -t {threads} 2>{log}"""

rule fanc_to_juicer:
    input:
        pairs = get_pairs_files,
        jar = "resources/juicer/juicer_tools.2.20.00.jar"
    output:
        "results/juicer/{sample_group}.{enzyme}.{fragments}.juicer.hic"
    params:
        files = lambda wc, input: " ".join(input.pairs) if isinstance(input.pairs,list) else input.pairs
    log:
        "logs/fanc/{sample_group}.{enzyme}.{fragments}.to_juicer.log"
    conda:
        "../envs/fanc.yaml"
    threads: 24
    shell:
        "fanc to-juicer {params.files} {output} --juicer-tools-jar {input.jar} 2>{log}"

rule fanc_to_cooler:
    input:
        hic = "results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic"
    output:
        "results/cooler/{sample_group}.{enzyme}.{fragments}-{resolution}.mcool"
    log:
        "logs/fanc/{sample_group}.{enzyme}.{fragments}-{resolution}.to_cooler.log"
    conda:
        "../envs/fanc.yaml"
    threads: 24
    shell:
        "fanc to-cooler {input} {output} -t {threads} 2>{log}"