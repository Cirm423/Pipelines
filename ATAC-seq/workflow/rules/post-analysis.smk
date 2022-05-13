rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        lambda wildcards, resources: f"-v {'' if config['single_end'] else '-pe -seg_len 100000'} -seed 1 {'-D' if resources.attempt > 1 else ''}" 
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
         bam="results/filtered/{sample}.sorted.bam",
         ref=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa"
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