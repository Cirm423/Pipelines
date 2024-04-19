rule plot_fingerprint:
    input:
        bam_files=["results/bamtools_filtered/{sample}.sorted.bam", "results/bamtools_filtered/{control}.sorted.bam"],
        bam_idx=["results/bamtools_filtered/{sample}.sorted.bam.bai", "results/bamtools_filtered/{control}.sorted.bam.bai"],
        jsd_sample="results/bamtools_filtered/{control}.sorted.bam",
        stats=expand("results/{step}/{{sample}}.sorted.{step}.stats.txt",
            step="bamtools_filtered")
    output:  #ToDo: add description to report caption
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint=report("results/deeptools/{sample}-{control}.plot_fingerprint.pdf", caption="../report/plot_fingerprint_deeptools.rst", category="QC"),
        counts="results/deeptools/{sample}-{control}.fingerprint_counts.txt",
        qc_metrics="results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--skipZeros ",
        "--numberOfSamples 500000 ", # ToDo: to config?
        lambda w, input:
            "{se_option}{fragment_size}".format(
                se_option="--extendReads " if config["single_end"] else "",
                # Estimated fragment size used to extend single-end reads
                fragment_size=
                "$(grep ^SN {stats} | "
                "cut -f 2- | "
                "grep -m1 'average length:' | "
                "awk '{{print $NF}}') ".format(
                    stats=input.stats)
                if config["single_end"] else ""
            )
    threads:
        8
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools/plotfingerprint.py"

rule prep_seacr:
    output:
        temp("SEACR_1.3.R")
    conda:
        "../envs/seacr.yaml"
    shell:
        "ln -s $(which SEACR_1.3.R) ./SEACR_1.3.R"

rule seacr_callpeak_stringent:
    input:
        "SEACR_1.3.R",
        sample=f"results/bed_graph/{{sample}}{'_normalized' if config['params']['callpeak']['spike'] else ''}.bedgraph",
        control=f"results/bed_graph/{{control}}{'_normalized' if config['params']['callpeak']['spike'] else ''}.bedgraph",
    output:
        "results/seacr_callpeak/{sample}-{control}.stringent.bed"
    params:
        extra=f"{'non' if config['params']['callpeak']['spike'] else 'norm'} stringent",
        out_prefix = lambda wc, output: output[0].split(".stringent.bed")[0]
    log:
        "logs/seacr/{sample}-{control}.log"
    conda:
        "../envs/seacr.yaml"
    shell:
        "bash SEACR_1.3.sh {input.sample} {input.control} {params.extra} {params.out_prefix} 2>{log}"

rule seacr_callpeak_relaxed:
    input:
        "SEACR_1.3.R",
        sample=f"results/bed_graph/{{sample}}{'_normalized' if config['params']['callpeak']['spike'] else ''}.bedgraph",
        control=f"results/bed_graph/{{control}}{'_normalized' if config['params']['callpeak']['spike'] else ''}.bedgraph",
    output:
        "results/seacr_callpeak/{sample}-{control}.relaxed.bed"
    params:
        extra=f"{'non' if config['params']['callpeak']['spike'] else 'norm'} relaxed",
        out_prefix = lambda wc, output: output[0].split(".relaxed.bed")[0]
    log:
        "logs/seacr/{sample}-{control}.log"
    conda:
        "../envs/seacr.yaml"
    shell:
        "bash SEACR_1.3.sh {input.sample} {input.control} {params.extra} {params.out_prefix} 2>{log}"

rule peaks_count:
    input:
        peaks=f"results/seacr_callpeak/{{sample}}-{{control}}.{config['params']['peak-analysis']}.bed"
    output:
        "results/seacr_callpeak/peaks_count/{sample}-{control}.peaks_count.tsv"
    log:
        "logs/seacr_callpeak/peaks_count/{sample}-{control}.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | "
        " wc -l | "
        " gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}_peaks\", $1}}' "
        " > {output} 2> {log}"

rule sm_report_peaks_count_plot:
    input:
        get_peaks_count_plot_input()
    output:
        report("results/seacr_callpeak/plots/plot_peaks_count.pdf", caption="../report/plot_peaks_count_macs2.rst", category="CallPeaks")
    log:
        "logs/seacr_callpeak/plot_peaks_count.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_peaks_count_macs2.R"

rule bedtools_intersect:
    input:
        left="results/bamtools_filtered/{sample}.sorted.bam",
        right=f"results/seacr_callpeak/{{sample}}-{{control}}.{config['params']['peak-analysis']}.bed"
    output:
        "results/bedtools_intersect/{sample}-{control}.intersected.bed"
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/bedtools/intersect/{sample}-{control}.intersected.log"
    wrapper:
        "v1.3.1/bio/bedtools/intersect"

rule frip_score:
    input:
        intersect="results/bedtools_intersect/{sample}-{control}.intersected.bed",
        flagstats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat", step= "bamtools_filtered")
    output:
        "results/bedtools_intersect/{sample}-{control}.peaks_frip.tsv"
    log:
        "logs/bedtools/intersect/{sample}-{control}.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep -m 1 'mapped (' {input.flagstats} | "
        " gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) "
        " -v OFS='\t' "
        " '{{print \"{wildcards.sample}-{wildcards.control}_peaks\", a/$1}}' "
        " > {output} 2> {log}"

rule sm_rep_frip_score:
    input:
        get_frip_score_input()
    output:
        report("results/seacr_callpeak/plots/plot_peaks_frip_score.pdf", caption="../report/plot_frip_score_macs2_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"
#May fail
rule create_igv_peaks:
    input:
        f"results/seacr_callpeak/{{sample}}-{{control}}.{config['params']['peak-analysis']}.bed"
    output:
        "results/IGV/seacr_callpeak/merged_library.{sample}-{control}.peaks.igv.txt"
    params:
        f"{config['params']['peak-analysis']}.bed"
    log:
        "logs/igv/create_igv_peaks/merged_library.{sample}-{control}.peaks.log"
    shell:
        " find {input} -type f -name '*.{params}' -exec echo -e 'results/IGV/seacr_callpeak/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"

rule homer_annotatepeaks:
    input:
        peaks=f"results/seacr_callpeak/{{sample}}-{{control}}.{config['params']['peak-analysis']}.bed",
        genome=f"{assembly_path}{assembly}.fa",
        gtf=f"{assembly_path}{assembly}.annotation.gtf"
    output:
        annotations="results/homer/annotate_peaks/{sample}-{control}.peaks.annotatePeaks.txt"
    threads:
        2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_peaks/{sample}-{control}.log"
    wrapper:
        "v1.3.1/bio/homer/annotatePeaks"

rule plot_homer_annotatepeaks:
    input:
        get_plot_homer_annotatepeaks_input()
    output:  #ToDo: add description to report caption
        summmary="results/homer/plots/plot_annotatepeaks_summary.txt",
        plot=report("results/homer/plots/plot_annotatepeaks.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "logs/homer/plot_annotatepeaks.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript workflow/scripts/plot_homer_annotatepeaks.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_sum_annotatepeaks:
    input:
        "results/homer/plots/plot_annotatepeaks_summary.txt"
    output:
        report("results/homer/plots/plot_annotatepeaks_summary.pdf", caption="../report/plot_annotatepeaks_summary_homer.rst", category="CallPeaks")
    log:
        "logs/homer/plot_annotatepeaks_summary.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_annotatepeaks_summary_homer.R"
