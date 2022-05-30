import os

rule sort_bed:
    input:
        "results/seacr_callpeak/{sam_contr_peak}.bed.peaks"
    output:
        "results/seacr_callpeak/{sam_contr_peak}.sorted.bed.peaks"
    log:
        "logs/seacr_callpeak/sort/{sam_contr_peak}.log"
    threads: 8
    conda:
        "../envs/bedsort.yaml"
    shell:
        "bedSort {input} {output} 2> {log}"

rule bedtools_merge_peaks:
    input:
        get_macs2_peaks_ab_sorted
    output:
        "results/bedtools/merged/{antibody}.consensus_peaks.txt"
    params:
        extra="-c 2,3,4,5,6,7,7 -o collapse,collapse,collapse,collapse,collapse,collapse,count_distinct"
    log:
        "logs/bedtools/merged/{antibody}.consensus_peaks.log"
    wrapper:
        "v1.3.1/bio/bedtools/merge"

# rule macs2_merged_expand:
#     input:
#         "results/bedtools/merged/{antibody}.consensus_peaks.txt"
#     output:
#         bool_txt="results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.txt",
#         bool_intersect="results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.intersect.txt"
#     params:
#         sample_control_peak=lambda wildcards: get_sample_control_peak_combinations_list_ab(wildcards.antibody),
#         narrow_param="--is_narrow_peak" if config["params"]["peak-analysis"] == "narrow" else "",
#         min_reps_consensus=config["params"]["min-reps-consensus"]
#     log:
#         "logs/macs2_merged_expand/{antibody}.consensus_peaks.boolean.log"
#     script:
#         "../scripts/macs2_merged_expand.py"

rule filter_consensus_peaks:
    input:
        "results/bedtools/merged/{antibody}.consensus_peaks.txt"
    output:
        "results/seacr_merged/{antibody}.consensus_peaks.filtered.txt"
    params:
        f"' \$10 >= {config['params']['min-reps-consensus']} {{print \$0}}'"
    log:
        "results/seacr_merged/{antibody}.consensus_peaks.filter.log"
    shell:
        "awk {params} {input} > {output} 2>{log}"

rule create_consensus_bed:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.txt"
    output:
        "results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.bed"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/macs2_merged_expand/{antibody}.consensus_peaks.boolean.bed.log"
    shell:
        "gawk -v FS='\t' -v OFS='\t' 'FNR  > 1 {{ print $1, $2, $3, $4 \"0\", \"+\"}}' {input} > {output} 2> {log}"

rule create_consensus_saf:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.txt"
    output:
        "results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.saf"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/macs2_merged_expand/{antibody}.consensus_peaks.boolean.bed.log"
    shell:
        "$(echo -e 'GeneID\tChr\tStart\tEnd\tStrand' > {output} && "
        " gawk -v FS='\t' -v OFS='\t' 'FNR > 1 {{ print $4, $1, $2, $3,  \" + \" }}' {input} >> {output}) "
        " 2> {log}"

rule plot_peak_consensus:
    input:
        "results/bedtools/merged/{antibody}.consensus_peaks.txt"
    output:
       report("results/seacr_merged/plots/{antibody}.consensus_peaks.pdf", caption="../report/plot_consensus_peak_intersect.rst", category="ConsensusPeak")
    params:
        lambda wildcards, output: os.path.dirname(output[0])
    conda:
        "../envs/consensus_plot.yaml"
    shell:
        "Rscript workflow/scripts/consensus_peaks.py --peaks {input} --outpath {params} 2> {log}"

rule create_consensus_igv:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.bed"
    output:
        "results/IGV/consensus/merged_library.{antibody}.consensus_peaks.igv.txt"
    log:
        "logs/igv/consensus/merged_library.{antibody}.consensus_peaks.igv.log"
    shell:
        "find {input} -type f -name '*.consensus_peaks.boolean.bed' -exec echo -e 'results/IGV/consensus/{wildcards.antibody}/\"{{}}\"\t0,0,0' \; > {output} 2> {log}"

rule homer_consensus_annotatepeaks:
    input:
        peaks="results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.bed",
        genome=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf"
    output:
        annotations="results/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.annotatePeaks.txt"
    threads:
        2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.annotatePeaks.log"
    wrapper:
        "v1.3.1/bio/homer/annotatePeaks"

rule trim_homer_consensus_annotatepeaks:
    input:
        "results/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.annotatePeaks.txt"
    output:
        temp("results/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.annotatePeaks.trimmed.txt")
    log:
        "logs/homer/annotate_consensus_peaks/trimmed/{antibody}.consensus_peaks.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cut -f2- {input} | gawk 'NR==1; NR > 1 {{print $0 | \"sort -T '.' -k1,1 -k2,2n\"}}' | cut -f6- > {output}"

rule merge_bool_and_annotatepeaks:
    input:
        trim="results/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.annotatePeaks.trimmed.txt",
        bool="results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.txt"
    output:
        "results/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.boolean.annotatePeaks.txt"
    log:
        "logs/homer/annotate_consensus_peaks/{antibody}.consensus_peaks.boolean.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "paste {input.bool} {input.trim} > {output}"

rule feature_counts:
    input:
        sam=lambda wc: expand(["results/filtered/{sample}.sorted.bam", "results/filtered/{control}.sorted.bam"],
            sample=get_samples_of_antibody(wc.antibody),
            control=get_controls_of_antibody(wc.antibody)),
        annotation="results/macs2_merged_expand/{antibody}.consensus_peaks.boolean.saf"
    output:
        multiext("results/feature_counts/{antibody}.consensus_peaks",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        2
    params:
        extra="-F SAF -O --fracOverlap 0.2{pe_param}".format(pe_param="" if config["single_end"] else " -p --donotsort")
    log:
        "logs/feature_counts/{antibody}.consensus_peaks.featureCounts.log"
    wrapper:
        "v1.1.0/bio/subread/featurecounts"

rule featurecounts_modified_colnames:
    input:
        featurecounts="results/feature_counts/{antibody}.consensus_peaks.featureCounts",
        bam=expand("results/filtered/{sample}.sorted.bam", sample=samples.index),
        samples_file=config["samples"]
    output:
        "results/feature_counts/{antibody}.consensus_peaks_modified.featureCounts"
    params:
        ""
    log:
        "logs/feature_counts/{antibody}.consensus_peaks_modified.featureCounts.log"
    script:
        "../scripts/col_mod_featurecounts.py"

rule featurecounts_deseq2:
    input:
        "results/feature_counts/{antibody}.consensus_peaks_modified.featureCounts"
    output:
        dds="results/deseq2/dss_rld/{antibody}.consensus_peaks.dds.rld.RData",
        plot_pca=report("results/deseq2/plots/{antibody}.consensus_peaks.pca_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_pca.rst", category = "DESeq2"),
        plot_heatmap=report("results/deseq2/plots/{antibody}.consensus_peaks.heatmap_plot.pdf",  #ToDo: add description to report caption
            caption = "../report/plot_deseq2_heatmap.rst", category = "DESeq2"),
        pca_data="results/deseq2/pca_vals/{antibody}.consensus_peaks.pca.vals.txt",
        dist_data="results/deseq2/dists/{antibody}.consensus_peaks.sample.dists.txt",
        size_factors_rdata="results/deseq2/sizeFactors/{antibody}.consensus_peaks.sizeFactors.RData",
        size_factors_res="results/deseq2/sizeFactors/{antibody}.consensus_peaks.sizeFactors.sizeFactor.txt",
        results="results/deseq2/results/{antibody}.consensus_peaks.deseq2_results.txt",
        FDR_1_perc_res="results/deseq2/FDR/{antibody}.consensus_peaks.deseq2.FDR_0.01.results.txt",
        FDR_5_perc_res="results/deseq2/FDR/{antibody}.consensus_peaks.deseq2.FDR_0.05.results.txt",
        FDR_1_perc_bed="results/deseq2/FDR/{antibody}.consensus_peaks.deseq2.FDR_0.01.results.bed",
        FDR_5_perc_bed="results/deseq2/FDR/{antibody}.consensus_peaks.deseq2.FDR_0.05.results.bed",
        plot_FDR_1_perc_MA=report("results/deseq2/plots/FDR/{antibody}.consensus_peaks_FDR_0.01_MA_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_1_perc_MA.rst", category = "DESeq2-FDR"),
        plot_FDR_5_perc_MA=report("results/deseq2/plots/FDR/{antibody}.consensus_peaks_FDR_0.05_MA_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_5_perc_MA.rst", category = "DESeq2-FDR"),
        plot_FDR_1_perc_volcano=report("results/deseq2/plots/FDR/{antibody}.consensus_peaks_FDR_0.01_volcano_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_1_perc_volcano.rst", category = "DESeq2-FDR"),
        plot_FDR_5_perc_volcano=report("results/deseq2/plots/FDR/{antibody}.consensus_peaks_FDR_0.05_volcano_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_5_perc_volcano.rst", category = "DESeq2-FDR"),
        plot_sample_corr_heatmap=report("results/deseq2/plots/{antibody}.consensus_peaks_sample_corr_heatmap.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_sample_corr_heatmap.rst", category = "DESeq2"),
        plot_scatter=report("results/deseq2/plots/{antibody}.consensus_peaks_scatter_plots.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_scatter.rst", category = "DESeq2")
    threads:
        2
    params:
        vst = config["params"]["deseq2"]["vst"]
    log:
        "logs/deseq2/{antibody}.consensus_peaks.featureCounts.log"
    conda:
        "../envs/featurecounts_deseq2.yaml"
    script:
        "../scripts/featurecounts_deseq2.R"

rule create_deseq2_igv:
    input:
        "results/deseq2/results/{antibody}.consensus_peaks.deseq2.FDR_0.05.results.bed"
    output:
        "results/IGV/consensus/merged_library.{antibody}.consensus_peaks.deseq2.FDR_0.05.igv.txt"
    log:
        "logs/igv/consensus/merged_library.{antibody}.consensus_peaks.deseq2.FDR_0.05.igv.log"
    shell:
        "find {input} -type f -name '*.consensus_{wildcards.peak}-peaks.deseq2.FDR_0.05.results.bed' -exec echo -e 'results/IGV/consensus/{wildcards.antibody}/deseq2/\"{{}}\"\t255,0,0' \; > {output} 2> {log}"
