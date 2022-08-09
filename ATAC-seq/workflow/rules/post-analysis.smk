rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        lambda wildcards, resources: f"-v {'' if config['single_end'] else '-pe -seg_len 1000000000'} -seed 1 {'-D' if resources.attempt > 1 else ''}" 
    log:
        "logs/preseq/{sample}.log"
    resources:
        attempt = lambda wildcards, attempt: attempt
    conda:
        "../envs/preseq.yaml"
    script:
        "../scripts/preseq.py"

rule collect_multiple_metrics:
    input:
         bam="results/bamtools_filtered/{sample}.sorted.bam",
         ref=f"{assembly_path}{assembly}.fa"
    output: #ToDo: add descriptions to report captions
        # Through the output file extensions the different tools for the metrics can be selected
        # so that it is not necessary to specify them under params with the "PROGRAM" option.
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html.
        [
           "{path}{sample}.insert_size_metrics",
            report(
                "{path}{sample}.insert_size_histogram.pdf",
                caption="../report/plot_insert_size_histogram_picard_mm.rst",
                category="Multiple Metrics (picard)"
            )
        ] if not config["single_end"] else [],
        multiext("{path}{sample}",
                 ".alignment_summary_metrics",
                 ".base_distribution_by_cycle_metrics",
                 ".quality_by_cycle_metrics",
                 ".quality_distribution_metrics",
                 ),
        report("{path}{sample}.base_distribution_by_cycle.pdf", caption="../report/plot_base_distribution_by_cycle_picard_mm.rst", category="Multiple Metrics (picard)"),
        report("{path}{sample}.quality_by_cycle.pdf", caption="../report/plot_quality_by_cycle_picard_mm.rst", category="Multiple Metrics (picard)"),
        report("{path}{sample}.quality_distribution.pdf", caption="../report/plot_quality_distribution_picard_mm.rst", category="Multiple Metrics (picard)")
    resources:
        # This parameter (default 3 GB) can be used to limit the total resources a pipeline is allowed to use, see:
        #     https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources
        mem_gb=3
    log:
        "logs/picard/{path}{sample}.log"
    params:
        # optional parameters, TODO: move to config.yaml and load from there
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
    conda:
        "../envs/picard.yaml"
    script:
        "../scripts/picard_metrics.py"

rule genomecov:
    input:
        "results/bamtools_filtered/{sample}.sorted.bam",
        flag_stats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat",
            step= "bamtools_filtered"),
        stats=expand("results/{step}/{{sample}}.sorted.{step}.stats.txt",
            step= "bamtools_filtered"),
    output:
        "results/bed_graph/{sample}.bedgraph"
    log:
        "logs/bed_graph/{sample}.log"
    params:
        lambda w, input:
            "-bg -scale $(grep -m 1 'mapped (' {flagstats_file} | awk '{{print 1000000/$1}}') {pe_fragment} {extend}".format(
            flagstats_file=input.flag_stats,
            pe_fragment="" if config["single_end"] else "-pc",
            # Estimated fragment size used to extend single-end reads
            extend=
                "-fs $(grep ^SN {stats} | "
                "cut -f 2- | "
                "grep -m1 'average length:' | "
                "awk '{{print $NF}}')".format(
                stats=input.stats)
            if config["single_end"] else ""
        )
    wrapper:
        "v1.3.1/bio/bedtools/genomecov"

rule sort_genomecov:
    input:
        "results/bed_graph/{sample}.bedgraph"
    output:
        "results/bed_graph/{sample}.sorted.bedgraph"
    log:
        "logs/sort_genomecov/{sample}.log"
    threads: 8
    conda:
        "../envs/bedsort.yaml"
    shell:
        "bedSort {input} {output} 2> {log}"

rule bedGraphToBigWig:
    input:
        bedGraph="results/bed_graph/{sample}.sorted.bedgraph",
        chromsizes=f"{assembly_path}{assembly}.chrom.sizes"
    output:
        "results/big_wig/{sample}.bigWig"
    log:
        "logs/big_wig/{sample}.log"
    params:
        ""
    wrapper:
        "v1.3.1/bio/ucsc/bedGraphToBigWig"

rule create_igv_bigwig:
    input:
        f"{assembly_path}{assembly}.annotation.bed",
        expand("results/big_wig/{sample}.bigWig", sample=samples.index)
    output:
        "results/IGV/big_wig/merged_library.bigWig.igv.txt"
    log:
        "logs/igv/create_igv_bigwig.log"
    shell:
        "find {input} -type f -name '*.bigWig' -exec echo -e 'results/IGV/big_wig/\"{{}}\"\t0,0,178' \;  > {output} 2> {log}"