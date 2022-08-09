#Genrich requires the bam files to be sorted by queryname instead of coordinate like all the stats and qc stuff, make temp files for that.
rule genrich_sort:
    input:
        "results/bamtools_filtered/{sample}.sorted.bam"
    output:
        temp("results/genrich/{sample}.sorted.bam")
    params:
        extra="-n"
    log:
        "logs/genrich/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"


#Genrich manages shifting, Removal of mitochondrial reads, Removal of PCR duplicates, Analysis of multimapping reads. The blacklisted regions are removed previously.
rule genrich:
    input:
        get_genrich_input
    output:
        peak = "results/genrich/{group}.narrowPeak",
        bed = "results/genrich/{group}.bed",
    log:
        "logs/genrich/{group}.log"
    params:
        samples = lambda wc: ",".join(get_samples_of_group(wc.group)),
        controls = lambda wc: "" if len(get_controls_of_group(wc.group)) == 0 else "-c " + ",".join(get_controls_of_group(wc.group)),
        #f Output bedgraph for visualizacion, r remove pcr duplicates, -e for excluding chromosomes, j is ATAC mode
        base = lambda wc, input, output: f"-f {output.bed} -r -j -v",
        #Tells genrich to consider unpaired alignments if single end
        single = "-y" if config["single_end"] else "",
        excl = "" if not config["params"]["callpeak"]["chromosome"] else f"-e {config['params']['callpeak']['chromosome']}",
        p_value = "-p {}".format(config["params"]["callpeak"]["p-value"]) if config["params"]["callpeak"]["p-value"] else "",
        q_value = "-q {}".format(config["params"]["callpeak"]["q-value"]) if config["params"]["callpeak"]["q-value"] else "",
    conda:
        "../envs/genrich.yaml"
    threads: 12
    shell:
        "Genrich {params.base} {params.single} {params.excl} {params.p_value} {params.q_value} -t {params.samples} {params.controls} -o {output.peak} 2>{log}"

 #The bedgraph output above has to be processed into 4 columns (first 3 + a datavalue of choice) for visualization  

#Only Frip score of treatment samples, ignoring controls or inputs.

rule peaks_count:
    input:
        peaks="results/genrich/{group}.narrowPeak"
    output:
        "results/genrich/peaks_count/{group}.narrow.peaks_count.tsv"
    log:
        "logs/genrich/peaks_count/{group}.narrow.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | "
        " wc -l | "
        " gawk -v OFS='\t' '{{print \"{wildcards.group}_narrow_peaks\", $1}}' "
        " > {output} 2> {log}"

rule sm_report_peaks_count_plot:
    input:
        expand("results/genrich/peaks_count/{group}.narrow.peaks_count.tsv", group = groups)
    output:
        report("results/genrich/plots/plot_narrow_peaks_count.pdf", caption="../report/plot_peaks_count_genrich.rst", category="CallPeaks")
    log:
        "logs/genrich/plot_narrow_peaks_count.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_peaks_count_genrich.R"

rule bedtools_intersect:
    input:
        left="results/bamtools_filtered/{sample}.sorted.bam",
        right= lambda wc: expand("results/genrich/{group}.narrowPeak", group = samples[samples.index == wc.sample]['group'])
    output:
        "results/bedtools_intersect/{sample}.intersected.bed"
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/bedtools/intersect/{sample}.intersected.log"
    wrapper:
        "v1.3.1/bio/bedtools/intersect"

rule frip_score:
    input:
        intersect="results/bedtools_intersect/{sample}.intersected.bed",
        flagstats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat", step= "bamtools_filtered")
    output:
        "results/bedtools_intersect/{sample}.narrow.peaks_frip.tsv"
    log:
        "logs/bedtools/intersect/{sample}.narrow.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep -m 1 'mapped (' {input.flagstats} | "
        " gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) "
        " -v OFS='\t' "
        " '{{print \"{wildcards.sample}_narrow_peaks\", a/$1}}' "
        " > {output} 2> {log}"

rule sm_rep_frip_score:
    input:
        expand("results/bedtools_intersect/{sample}.narrow.peaks_frip.tsv",sample = samples.loc[samples['group'].isin(groups)].index)
    output:
        report("results/genrich/plots/plot_narrow_peaks_frip_score.pdf", caption="../report/plot_frip_score_genrich_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_narrow_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"

rule create_igv_peaks:
    input:
        "results/genrich/{group}.narrowPeak"
    output:
        "results/IGV/genrich_peaks/merged_library.{group}.narrow_peaks.igv.txt"
    log:
        "logs/igv/genrich_peaks/merged_library.{group}.narrow_peaks.log"
    shell:
        " find {input} -type f -name '*.narrowPeak' -exec echo -e 'results/IGV/genrich_peaks/narrow/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"

rule homer_annotatepeaks:
    input:
        peaks="results/genrich/{group}.narrowPeak",
        genome=f"{assembly_path}{assembly}.fa",
        gtf=f"{assembly_path}{assembly}.annotation.gtf"
    output:
        annotations="results/homer/annotate_peaks/{group}.narrow_peaks.annotatePeaks.txt"
    threads:
        2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_peaks/{group}.narrow.log"
    wrapper:
        "v1.3.1/bio/homer/annotatePeaks"

rule plot_macs_qc:
    input:
        expand("results/genrich/{group}.narrowPeak", group = groups)
    output:  #ToDo: add description to report caption
        summmary="results/genrich/plots/plot_narrow_peaks_genrich_summary.txt",
        plot=report("results/genrich/plots/plot_narrow_peaks_genrich.pdf", caption="../report/plot_genrich_qc.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        groups = ','.join(groups)
    log:
        "logs/genrich/plot_narrow_peaks_genrich.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript workflow/scripts/plot_macs_qc.R -i {params.input} -s {params.groups}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_homer_annotatepeaks:
    input:
        expand("results/homer/annotate_peaks/{group}.narrow_peaks.annotatePeaks.txt", group = groups)
    output:  #ToDo: add description to report caption
        summmary="results/homer/plots/plot_narrow_annotatepeaks_summary.txt",
        plot=report("results/homer/plots/plot_narrow_annotatepeaks.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        groups = ','.join(groups)
    log:
        "logs/homer/plot_narrow_annotatepeaks.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript workflow/scripts/plot_homer_annotatepeaks.R -i {params.input} -s {params.groups}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_sum_annotatepeaks:
    input:
        "results/homer/plots/plot_narrow_annotatepeaks_summary.txt"
    output:
        report("results/homer/plots/plot_narrow_annotatepeaks_summary.pdf", caption="../report/plot_annotatepeaks_summary_homer.rst", category="CallPeaks")
    log:
        "logs/homer/plot_narrow_annotatepeaks_summary.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_annotatepeaks_summary_homer.R"
