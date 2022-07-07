rule preseq_lc_extrap:
    input:
        get_dedup_bam
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        lambda wildcards, resources: f"-v {'' if config['single_end'] else '-pe -seg_len 1000000000'} -seed 1 {'-D' if resources.attempt > 1 else ''}" 
    log:
        "logs/preseq/{sample}.log"
    resources:
        attempt = lambda wildcards, attempt: attempt
    threads: 8
    conda:
        "../envs/preseq.yaml"
    script:
        "../scripts/preseq.py"

rule bismark2report_pe:
    input:
        alignment_report="results/bismark_mapped/{sample}_PE_report.txt",
        dedup_report="results/bismark_mapped/{sample}_pe.deduplication_report.txt",
        mbias_report="results/bismark/meth/{sample}-pe.M-bias.txt",
        splitting_report="results/bismark/meth/{sample}_pe.deduplicated_splitting_report.txt"
    output:
        html="results/qc/bismark/{sample}_pe.bismark2report.html",
    log:
        "logs/qc/bismark/{sample}_pe.bismark2report.html.log",
    params:
        skip_optional_reports=True
    wrapper:
        "v1.7.0/bio/bismark/bismark2report"

rule bismark2report_se:
    input:
        alignment_report="results/bismark_mapped/{sample}_SE_report.txt",
        dedup_report="results/bismark_mapped/{sample}.deduplication_report.txt",
        mbias_report="results/bismark/meth/{sample}-se.M-bias.txt",
        splitting_report="results/bismark/meth/{sample}.deduplicated_splitting_report.txt"
    output:
        html="results/qc/bismark/{sample}_se.bismark2report.html",
    log:
        "logs/qc/bismark/{sample}_se.bismark2report.html.log",
    params:
        skip_optional_reports=True
    wrapper:
        "v1.7.0/bio/bismark/bismark2report"

rule bismark2summary_prepare_symlinks:
    input:
        get_sample_splitting_reports,
    output:
        temp(f"results/bismark_mapped/{{sample}}{'' if config['single_end'] else '_pe'}.deduplicated_splitting_report.txt"),
    log:
        "qc/bismark/{sample}_prepare_symlinks.symlinks.log"
    run:
        wd = os.getcwd()
        shell("echo 'Making symlinks' > {log}")
        for source, target in zip(input, output):
           target_dir = os.path.dirname(target)
           target_name = os.path.basename(target)
           log_path = os.path.join(wd, log[0])
           abs_src_path = os.path.abspath(source)
           shell("cd {target_dir} && ln -f -s {abs_src_path} {target_name} >> {log_path} 2>&1")

        shell("echo 'Done' >> {log}")

rule bismark2summary:
    input:
        bam=get_bismark_bams,

        # Bismark `bismark2summary` discovers reports automatically based
        # on files available in bam file containing folder
        #
        # If your per BAM file reports aren't in the same folder
        # you will need an additional task which symlinks all reports
        # (E.g. your splitting report generated by `bismark_methylation_extractor`
        # tool is in `meth` folder, and alignment related reports in `bams` folder)

        # These dependencies are here just to ensure that corresponding rules
        # has already finished at rule execution time, otherwise some reports
        # will be missing.
        dependencies=get_sample_splitting_reports_ln,
    output:
        html="qc/bismark/bismark2summary.html",
        txt="qc/bismark/bismark2summary.txt"
    log:
        "logs/qc/bismark/bismark2summary.log"
    wrapper:
        "v1.7.0/bio/bismark/bismark2summary"
