rule preseq_lc_extrap:
    input:
        "results/picard_dedup/{sample}.sorted.bam"
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

#Old rule for ChIP-seq, here using the SEACR bedgraphs instead
# rule genomecov:
#     input:
#         "results/bamtools_filtered/{sample}.sorted.bam",
#         flag_stats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat",
#             step= "bamtools_filtered" if config["single_end"]
#             else "orph_rm_pe"),
#         stats=expand("results/{step}/{{sample}}.sorted.{step}.stats.txt",
#             step= "bamtools_filtered" if config["single_end"]
#             else "orph_rm_pe"),
#     output:
#         "results/bed_graph/{sample}.bedgraph"
#     log:
#         "logs/bed_graph/{sample}.log"
#     params:
#         lambda w, input:
#             "-bg -scale $(grep -m 1 'mapped (' {flagstats_file} | awk '{{print 1000000/$1}}') {pe_fragment} {extend}".format(
#             flagstats_file=input.flag_stats,
#             pe_fragment="" if config["single_end"] else "-pc",
#             # Estimated fragment size used to extend single-end reads
#             extend=
#                 "-fs $(grep ^SN {stats} | "
#                 "cut -f 2- | "
#                 "grep -m1 'average length:' | "
#                 "awk '{{print $NF}}')".format(
#                 stats=input.stats)
#             if config["single_end"] else ""
#         )
#     wrapper:
#         "v1.3.1/bio/bedtools/genomecov"

rule sort_genomecov:
    input:
        "results/bed_graph/{sample}_normalized.bedgraph"
    output:
        "results/bed_graph/{sample}_normalized.sorted.bedgraph"
    log:
        "logs/sort_genomecov/{sample}.log"
    threads: 8
    conda:
        "../envs/bedsort.yaml"
    shell:
        "bedSort {input} {output} 2> {log}"

rule bedGraphToBigWig:
    input:
        bedGraph="results/bed_graph/{sample}_normalized.sorted.bedgraph",
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

rule create_region:
    input:
        get_macs2_peaks()
    output:
        "results/seacr_callpeak/Peak_region.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output}"

rule compute_matrix:
    input:
         #bed=f"{assembly_path}{assembly}.annotation.bed",
         bed="results/seacr_callpeak/Peak_region.bed",
         bigwig=expand("results/big_wig/{sample}.bigWig", sample=samples.index)
    output:
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz="results/deeptools/matrix_files/matrix.gz",
        matrix_tab="results/deeptools/matrix_files/matrix.tab"
    log:
        "logs/deeptools/compute_matrix.log"
    threads: 24
    params:
        command="scale-regions",
        extra="--numberOfProcessors 24 " #Match this to threads above
              "--regionBodyLength 1000 "
              "--beforeRegionStartLength 3000 "
              "--afterRegionStartLength 3000 "
              "--missingDataAsZero " # added to prevent black output in the heatmap (plot_heatmap rule) https://github.com/deeptools/deepTools/issues/793
              "--skipZeros "
              "--smartLabels "       
    conda: 
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools/computeMatrix.py"

rule plot_profile:
    input:
         "results/deeptools/matrix_files/matrix.gz"
    output: #ToDo: add description to report caption
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotprofile.html.
        plot_img=report("results/deeptools/plot_profile.pdf", caption="../report/plot_profile_deeptools.rst", category="GenomicRegions"),
        data="results/deeptools/plot_profile_data.tab"
    log:
        "logs/deeptools/plot_profile.log"
    params:
        ""
    threads: 24
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools/plotProfile.py"

rule plot_heatmap:
    input:
         "results/deeptools/matrix_files/matrix.gz"
    output:  #ToDo: add description to report caption
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotheatmap.html.
        heatmap_img=report("results/deeptools/heatmap.pdf", caption="../report/plot_heatmap_deeptools.rst", category="Heatmaps"),
        heatmap_matrix="results/deeptools/heatmap_matrix.tab"
    log:
        "logs/deeptools/heatmap.log"
    params:
        ""
    threads: 24
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools/plotheatmap.py"

rule phantompeakqualtools:
    input:
        "results/bamtools_filtered/{sample}.sorted.bam"
    output:  #ToDo: add description to report caption
        res_phantom="results/phantompeakqualtools/{sample}.phantompeak.spp.out",
        r_data="results/phantompeakqualtools/{sample}.phantompeak.Rdata",
        plot=report("results/phantompeakqualtools/{sample}.phantompeak.pdf", caption="../report/plot_phantompeak_phantompeakqualtools.rst", category="Phantompeak")
    threads:
        8
    log:
        "logs/phantompeakqualtools/{sample}.phantompeak.log"
    conda:
        "../envs/phantompeakqualtools.yaml"
    shell:
        "( Rscript -e \"library(caTools); source('workflow/scripts/run_spp.R')\" "
        "  -c={input} -savp={output.plot} -savd={output.r_data} "
        "  -out={output.res_phantom} -p={threads} 2>&1 ) >{log}"

rule phantompeak_correlation:
    input:
        data="results/phantompeakqualtools/{sample}.phantompeak.Rdata",
        header="workflow/header/spp_corr_header.txt"
    output:
        "results/phantompeakqualtools/{sample}.spp_correlation_mqc.tsv"
    log:
        "logs/phantompeakqualtools/correlation/{sample}.spp_corr.log"
    conda:
        "../envs/phantom_corr.yaml"
    script:
        "../scripts/phantompeak_correlation.R"

rule phantompeak_multiqc:
    # NSC (Normalized strand cross-correlation) and RSC (relative strand cross-correlation) metrics use cross-correlation
    input:
        data="results/phantompeakqualtools/{sample}.phantompeak.spp.out",
        nsc_header="workflow/header/nsc_header.txt",
        rsc_header="workflow/header/rsc_header.txt"
    output:
        nsc="results/phantompeakqualtools/{sample}.spp_nsc_mqc.tsv",
        rsc="results/phantompeakqualtools/{sample}.spp_rsc_mqc.tsv"
    log:
        "logs/phantompeakqualtools/correlation/{sample}.nsc_rsc.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "( gawk -v OFS='\t' '{{print $1, $9}}' {input.data} | cat {input.nsc_header} - > {output.nsc} && "
        "  gawk -v OFS='\t' '{{print $1, $10}}' {input.data} | cat {input.rsc_header} - > {output.rsc} 2>&1 ) >{log}"
