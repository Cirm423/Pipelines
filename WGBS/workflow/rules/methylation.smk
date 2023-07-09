rule methyldackel:
    input:
        bam = "results/picard_dedup/{sample}.bam",
        bai = "results/picard_dedup/{sample}.bam.bai",
        fa = f"{assembly_path}{assembly}.fa",
        fai = f"{assembly_path}{assembly}.fa.fai",
    output:
        expand("results/methyldackel/{{sample}}_{meth}{suffix}",
            meth = ["CpG","CHG","CHH"] if config["params"]["methyldackel"]["comprehensive"] else ["CpG"],
            suffix = ".methylKit" if config["params"]["methyldackel"]["methyl_kit"] else ".bedGraph"
        ),
        txt = "results/methyldackel/{sample}_methyldackel.txt",
        svg_OT = report("results/methyldackel/{sample}_meth_bias_OT.svg", category="Methylation"),
        svg_OB = report("results/methyldackel/{sample}_meth_bias_OB.svg", category="Methylation")
    params:
        prefix = "results/methyldackel/{sample}",
        prefix_mbias = "results/methyldackel/{sample}_meth_bias",
        comprehensive = "--CHG --CHH" if config["params"]["methyldackel"]["comprehensive"] else "",
        min_depth = config["params"]["methyldackel"]["min_depth"] if config["params"]["methyldackel"]["min_depth"] > 0 else "",
        ignore_flags = "--ignoreFlags" if config["params"]["methyldackel"]["ignore_flags"] else "",
        methyl_kit = "--methylKit" if config["params"]["methyldackel"]["methyl_kit"] else "",
        extra_extract = config["params"]["methyldackel"]["extra_extract"],
        extra_mbias = config["params"]["methyldackel"]["extra_mbias"],
    log:
        "logs/methyldackel/{sample}_methylation.log"
    conda:
        "../envs/methyldackel.yaml"
    threads: 24
    shell:
        "MethylDackel extract -@ {threads} {params.comprehensive} {params.ignore_flags} {params.methyl_kit} {params.min_depth} {input.fa} {input.bam} {params.extra_extract} -o {params.prefix} 2>{log} && "
        "MethylDackel mbias -@ {threads} {params.comprehensive} {params.ignore_flags} {params.extra_mbias} {input.fa} {input.bam} {params.prefix_mbias} --txt > {output.txt} 2>>{log}"

rule bismark_methylation_extractor_pe:
    input: "results/bismark_mapped/{sample}_pe.deduplicated.bam"
    output:
        mbias_r1="results/qc/bismark/{sample}-pe.M-bias_R1.png",
        # Only for PE BAMS:
        mbias_r2="results/qc/bismark/{sample}-pe.M-bias_R2.png",

        mbias_report="results/bismark/meth/{sample}-pe.M-bias.txt",
        splitting_report=report("results/bismark/meth/{sample}_pe.deduplicated_splitting_report.txt", category = "Methylation"),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylome_CpG_cov="results/bismark/meth_cpg/{sample}-pe.bismark.cov.gz",
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph="results/bismark/meth_cpg/{sample}-pe.bedGraph.gz",

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="results/bismark/meth/CpG_context_{sample}-pe.txt.gz",
        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="results/bismark/meth/CHG_context_{sample}-pe.txt.gz",
        read_base_meth_state_chh="results/bismark/meth/CHH_context_{sample}-pe.txt.gz"
    log:
        "logs/bismark/{sample}-pe_methylaction_extraction.log"
    params:
        #These 2 params bellow may be used with other protocols
        # ignore_r2=2,
        # ignore_3prime_r2=2,
        output_dir="results/bismark/meth",  # optional output dir
        # optional params string, 8 threads only because bismark uses 3 * core processes
        extra=f"""--paired-end --no_overlap --gzip --multicore 8 --bedGraph --counts \\
        '--comprehensive ' \\
        --cutoff {config['params']['bismark']['extract']['cutoff']} \\
        {'--cytosine_report --genome_folder {} '.format(assembly_path) if config['params']['bismark']['extract']['cytosine_report'] else ''} \\
        {config['params']['bismark']['extract']['extra']}"""
    threads: 24
    conda:
        "../envs/bismark_meth.yaml"
    script:
        "../scripts/bismark_meth.py"

rule bismark_methylation_extractor_se:
    input: "results/bismark_mapped/{sample}_se.deduplicated.bam"
    output:
        mbias_r1="results/qc/bismark/{sample}-se.M-bias_R1.png",
        # Only for PE BAMS:
        # mbias_r2="qc/meth/{sample}.M-bias_R2.png",

        mbias_report="results/bismark/meth/{sample}-se.M-bias.txt",
        splitting_report=report("results/bismark/meth/{sample}.deduplicated_splitting_report.txt", category = "Methylation"),

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylome_CpG_cov="results/bismark/meth_cpg/{sample}-se.bismark.cov.gz",
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph="results/bismark/meth_cpg/{sample}-se.bedGraph.gz",

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="results/bismark/meth/CpG_context_{sample}-se.txt.gz" if config["params"]["bismark"]["extract"]["comprehensive"] else "",
        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="results/bismark/meth/CHG_context_{sample}-se.txt.gz" if config["params"]["bismark"]["extract"]["comprehensive"] else "",
        read_base_meth_state_chh="results/bismark/meth/CHH_context_{sample}-se.txt.gz" if config["params"]["bismark"]["extract"]["comprehensive"] else ""
    log:
        "logs/bismark/{sample}-se_methylaction_extraction.log"
    params:
        output_dir="results/bismark/meth",  # optional output dir
        # optional params string, 8 threads only because bismark uses 3 * core processes
        extra=f"""--single-end --gzip --multicore 8 --bedGraph --counts \\
        {' --comprehensive ' if config['params']['bismark']['extract']['comprehensive'] else ''} \\
        --cutoff {config['params']['bismark']['extract']['cutoff']} \\
        {' --cytosine_report --genome_folder {} '.format(assembly_path) if config['params']['bismark']['extract']['cytosine_report'] else ''} \\
        {config['params']['bismark']['extract']['extra']}"""
    threads: 24
    conda:
        "../envs/bismark_meth.yaml"
    script:
        "../scripts/bismark_meth.py"
