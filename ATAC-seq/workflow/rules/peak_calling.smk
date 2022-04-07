#Genrich manages shifting, Removal of mitochondrial reads, Removal of PCR duplicates, Analysis of multimapping reads.

rule genrich:
    input:
        samples = lambda wc: expand(["results/merged/{sample}.bam"],
            sample = get_samples_of_group(wc.group)),
        controls = lambda wc: expand(["results/merged/{control}.bam"],
            control = get_controls_of_group(wc.group)),
        blacklist = f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist_formated.sorted",
    output:
        peak = "results/genrich/{group}.narrowPeak",
        bed = "results/genrich/{group}.bed"
    log:
        "logs/genrich/{group}.log"
    params:
        #f Output bedgraph for visualizacion, r remove pcr duplicates, -e for excluding chromosomes, j is ATAC mode
        base = "-f -r -E {input.blacklist} -j -v",
        #Tells genrich to consider unpaired alignments if single end
        single = "-y" if config["single_end"] else "",
        excl = "" if not config["params"]["callpeak"]["chromosome"] else f"-e {config['params']['callpeak']['chromosome']}",
        p_value = "-p {}".format(config["params"]["callpeak"]["p-value"]) if config["params"]["callpeak"]["p-value"] else "",
        q_value = "-q {}".format(config["params"]["callpeak"]["q-value"]) if config["params"]["callpeak"]["q-value"] else "",
    conda:
        "../envs/genrich.yaml"
    shell:
        "Genrich -t {input.samples} -c {input.controls} -o {output.peak} {params.base} {params.single} {params.excl} {params.p_value} {params.q_value} 2>{log}"

 #The bedgraph output above has to be processed into 4 columns (first 3 + a datavalue of choice) for visualization  

rule merge_bams_group:
    input:
        lambda wc: expand("results/merged/{sample}.bam",
            sample = get_samples_of_group(wc.group)
        )
    output:
        temp("results/merged_group/{group}.bam")
    log:
        "logs/picard/mergebamfiles/{group}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate",
    wrapper:
        "v0.87.0/bio/picard/mergesamfiles"

rule samtools_sort:
    input:
        "results/merged_group/{group}.bam"
    output:
        temp("results/merged_group/{group}.sorted.bam")
    params:
        extra="-n"
    log:
        "logs/bamtools_filtered/{group}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule bedtools_intersect:
    input:
        left="results/merged_group/{group}.sorted.bam",
        right="results/genrich/{group}.narrowPeak"
    output:
        "results/bedtools_intersect/{group}.intersected.bed"
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/bedtools/intersect/{group}.intersected.log"
    wrapper:
        "v1.3.1/bio/bedtools/intersect"

rule frip_score:
    input:
        intersect="results/bedtools_intersect/{group}.intersected.bed",
        flagstats="results/merged_group/{group}.sorted.merged_group.flagstat"
    output:
        "results/bedtools_intersect/{group}.narrow.peaks_frip.tsv"
    log:
        "logs/bedtools/intersect/{group}.narrow.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep -m 1 'mapped (' {input.flagstats} | "
        " gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) "
        " -v OFS='\t' "
        " '{{print \"{wildcards.group}_narrow_peaks\", a/$1}}' "
        " > {output} 2> {log}"

rule sm_rep_frip_score:
    input:
        expand("results/bedtools_intersect/{group}.intersected.bed",group = groups)
    output:
        report("results/genrich/plots/plot_narrow_peaks_frip_score.pdf", caption="../report/plot_frip_score_genrich_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_narrow_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"

