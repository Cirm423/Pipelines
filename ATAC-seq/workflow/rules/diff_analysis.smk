rule bedtools_merge_narrow:
    input:
        expand("results/genrich/{group}.narrowPeak", group = groups)
    output:
        "results/bedtools/merged/consensus_narrow-peaks.txt"
    params:
        extra="-c {} -o {}".format( ','.join(map(str, list( range(2,11) ) ) ),
                                       ','.join( ["collapse"] * 9))
    log:
        "logs/bedtools/merged/consensus_peaks.log"
    wrapper:
        "v1.3.1/bio/bedtools/merge"

rule genrich_merged_expand:
    input:
        "results/bedtools/merged/consensus_narrow-peaks.txt"
    output:
        bool_txt="results/genrich_merged_expand/consensus_narrow-peaks.boolean.txt",
        bool_intersect="results/genrich_merged_expand/consensus_narrow-peaks.boolean.intersect.txt"
    params:
        sample_control_peak=expand("{group}.narrow", group = groups),
        narrow_param="--is_narrow_peak",
        min_reps_consensus=config["params"]["min-reps-consensus"]
    log:
        "logs/genrich_merged_expand/consensus_narrow-peaks.boolean.log"
    script:
        "../scripts/macs2_merged_expand.py"

rule create_consensus_bed:
    input:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.txt"
    output:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.bed"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/genrich_merged_expand/consensus_narrow-peaks.boolean.bed.log"
    shell:
        "gawk -v FS='\t' -v OFS='\t' 'FNR  > 1 {{ print $1, $2, $3, $4 \"0\", \"+\"}}' {input} > {output} 2> {log}"

rule create_consensus_saf:
    input:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.txt"
    output:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.saf"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/genrich_merged_expand/consensus_narrow-peaks.boolean.bed.log"
    shell:
        "$(echo -e 'GeneID\tChr\tStart\tEnd\tStrand' > {output} && "
        " gawk -v FS='\t' -v OFS='\t' 'FNR > 1 {{ print $4, $1, $2, $3,  \" + \" }}' {input} >> {output}) "
        " 2> {log}"

rule plot_peak_intersect:
    input:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.intersect.txt"
    output:
       report("results/genrich_merged_expand/plots/consensus_narrow-peaks.boolean.intersect.plot.pdf", caption="../report/plot_consensus_peak_intersect.rst", category="ConsensusPeak")
    conda:
        "../envs/consensus_plot.yaml"
    log:
        "logs/genrich_merged_expand/plots/consensus_narrow-peaks.boolean.intersect.plot.log"
    shell:
        "Rscript workflow/scripts/plot_peak_intersect.R -i {input} -o {output} 2> {log}"

rule create_consensus_igv:
    input:
        "results/genrich_merged_expand/consensus_narrow-peaks.boolean.bed"
    output:
        "results/IGV/consensus/merged_library.consensus_narrow-peaks.igv.txt"
    log:
        "logs/igv/consensus/merged_library.consensus_narrow-peaks.igv.log"
    shell:
        "find {input} -type f -name '*.consensus_narrow-peaks.boolean.bed' -exec echo -e 'results/IGV/consensus/\"{{}}\"\t0,0,0' \; > {output} 2> {log}"

rule homer_consensus_annotatepeaks:
    input:
        peaks="results/genrich_merged_expand/consensus_narrow-peaks.boolean.bed",
        genome=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf"
    output:
        annotations="results/homer/annotate_consensus_peaks/consensus_narrow-peaks.annotatePeaks.txt"
    threads:
        2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_consensus_peaks/consensus_narrow-peaks.annotatePeaks.log"
    wrapper:
        "v1.3.1/bio/homer/annotatePeaks"

rule trim_homer_consensus_annotatepeaks:
    input:
        "results/homer/annotate_consensus_peaks/consensus_narrow-peaks.annotatePeaks.txt"
    output:
        temp("results/homer/annotate_consensus_peaks/consensus_narrow-peaks.annotatePeaks.trimmed.txt")
    log:
        "logs/homer/annotate_consensus_peaks/trimmed/consensus_narrow-peaks.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cut -f2- {input} | gawk 'NR==1; NR > 1 {{print $0 | \"sort -T '.' -k1,1 -k2,2n\"}}' | cut -f6- > {output}"

rule merge_bool_and_annotatepeaks:
    input:
        trim="results/homer/annotate_consensus_peaks/consensus_narrow-peaks.annotatePeaks.trimmed.txt",
        bool="results/genrich_merged_expand/consensus_narrow-peaks.boolean.txt"
    output:
        "results/homer/annotate_consensus_peaks/consensus_narrow-peaks.boolean.annotatePeaks.txt"
    log:
        "logs/homer/annotate_consensus_peaks/consensus_narrow-peaks.boolean.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "paste {input.bool} {input.trim} > {output}"

rule feature_counts:
    input:
        sam=expand("results/filtered/{sample}.sorted.bam",
            sample = samples.index),
        annotation="results/genrich_merged_expand/consensus_narrow-peaks.boolean.saf"
    output:
        multiext("results/feature_counts/consensus_narrow-peaks",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        2
    params:
        extra="-F SAF -O --fracOverlap 0.2{pe_param}".format(pe_param="" if config["single_end"] else " -p --donotsort")
    log:
        "logs/feature_counts/consensus_narrow-peaks.featureCounts.log"
    wrapper:
        "v1.1.0/bio/subread/featurecounts"

rule featurecounts_modified_colnames:
    input:
        featurecounts="results/feature_counts/consensus_narrow-peaks.featureCounts",
        bam=expand("results/filtered/{sample}.sorted.bam", sample=samples.index),
        samples_file=config["samples"]
    output:
        "results/feature_counts/consensus_narrow-peaks_modified.featureCounts"
    params:
        ""
    log:
        "logs/feature_counts/consensus_narrow-peaks_modified.featureCounts.log"
    script:
        "../scripts/col_mod_featurecounts.py"

rule featurecounts_deseq2:
    input:
        "results/feature_counts/consensus_narrow-peaks_modified.featureCounts"
    output:
        dds="results/deseq2/dss_rld/consensus_narrow-peaks.dds.rld.RData",
        plot_pca=report("results/deseq2/plots/consensus_narrow-peaks.pca_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_pca.rst", category = "DESeq2"),
        plot_heatmap=report("results/deseq2/plots/consensus_narrow-peaks.heatmap_plot.pdf",  #ToDo: add description to report caption
            caption = "../report/plot_deseq2_heatmap.rst", category = "DESeq2"),
        pca_data="results/deseq2/pca_vals/consensus_narrow-peaks.pca.vals.txt",
        dist_data="results/deseq2/dists/consensus_narrow-peaks.sample.dists.txt",
        size_factors_rdata="results/deseq2/sizeFactors/consensus_narrow-peaks.sizeFactors.RData",
        size_factors_res="results/deseq2/sizeFactors/consensus_narrow-peaks.sizeFactors.sizeFactor.txt",
        results="results/deseq2/results/consensus_narrow-peaks.deseq2_results.txt",
        FDR_1_perc_res="results/deseq2/FDR/consensus_narrow-peaks.deseq2.FDR_0.01.results.txt",
        FDR_5_perc_res="results/deseq2/FDR/consensus_narrow-peaks.deseq2.FDR_0.05.results.txt",
        FDR_1_perc_bed="results/deseq2/FDR/consensus_narrow-peaks.deseq2.FDR_0.01.results.bed",
        FDR_5_perc_bed="results/deseq2/FDR/consensus_narrow-peaks.deseq2.FDR_0.05.results.bed",
        plot_FDR_1_perc_MA=report("results/deseq2/plots/FDR/consensus_narrow-peaks_FDR_0.01_MA_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_1_perc_MA.rst", category = "DESeq2-FDR"),
        plot_FDR_5_perc_MA=report("results/deseq2/plots/FDR/consensus_narrow-peaks_FDR_0.05_MA_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_5_perc_MA.rst", category = "DESeq2-FDR"),
        plot_FDR_1_perc_volcano=report("results/deseq2/plots/FDR/consensus_narrow-peaks_FDR_0.01_volcano_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_1_perc_volcano.rst", category = "DESeq2-FDR"),
        plot_FDR_5_perc_volcano=report("results/deseq2/plots/FDR/consensus_narrow-peaks_FDR_0.05_volcano_plot.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_FDR_5_perc_volcano.rst", category = "DESeq2-FDR"),
        plot_sample_corr_heatmap=report("results/deseq2/plots/consensus_narrow-peaks_sample_corr_heatmap.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_sample_corr_heatmap.rst", category = "DESeq2"),
        plot_scatter=report("results/deseq2/plots/consensus_narrow-peaks_scatter_plots.pdf", #ToDo: add description to report caption
            caption = "../report/plot_deseq2_scatter.rst", category = "DESeq2")
    threads:
        2
    params:
        vst = config["params"]["deseq2"]["vst"]
    log:
        "logs/deseq2/consensus_narrow-peaks.featureCounts.log"
    conda:
        "../envs/featurecounts_deseq2.yaml"
    script:
        "../scripts/featurecounts_deseq2.R"

rule create_deseq2_igv:
    input:
        "results/deseq2/results/consensus_narrow-peaks.deseq2.FDR_0.05.results.bed"
    output:
        "results/IGV/consensus/merged_library.consensus_narrow-peaks.deseq2.FDR_0.05.igv.txt"
    log:
        "logs/igv/consensus/merged_library.consensus_narrow-peaks.deseq2.FDR_0.05.igv.log"
    shell:
        "find {input} -type f -name '*.consensus_narrow-peaks.deseq2.FDR_0.05.results.bed' -exec echo -e 'results/IGV/consensus/deseq2/\"{{}}\"\t255,0,0' \; > {output} 2> {log}"