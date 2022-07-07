rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    threads: 6
    wrapper:
        "v1.3.1/bio/fastqc"

rule fastqc_trimmed:
    input:
        get_individual_trimmed_fastq
    output:
        html="results/qc/fastqc/trimmed_{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/trimmed_{sample}.{unit}.{read}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/trimmed_{sample}.{unit}.{read}.log"
    threads: 6
    wrapper:
        "v1.3.1/bio/fastqc"

rule sort_deduplicated:
    input:
        get_dedup_bam
    output:
        temp("results/qualimap/{sample}.sorted.bam")
    params:
        extra=""
    log:
        "logs/qualimap/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

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

rule qualimap:
    input:
        "results/qualimap/{sample}.sorted.bam"
    output:
        directory("results/qualimap/{sample}_qualimap")
    log:
        "logs/qualimap/{sample}.log"
    params:
        gcref = "-gd HUMAN" if "h" in assembly else "-gd MOUSE" if "m" in assembly else ""
    threads: 8
    conda:
        "../envs/qualimap.yaml"
    shell:
        "qualimap bamqc {params.gcref} -bam {input} -outdir {output} --collect-overlap-pairs -nt {threads} --java-mem-size=16G"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "results/qc/multiqc/multiqc.html"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/multiqc.py"