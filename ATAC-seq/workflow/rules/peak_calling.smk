#Genrich manages shifting, Removal of mitochondrial reads, Removal of PCR duplicates, Analysis of multimapping reads. The blacklisted regions are removed previously.

rule genrich:
    input:
        samples = lambda wc: expand(["results/filtered/{sample}.sorted.bam"],
            sample = get_samples_of_group(wc.group)),
        controls = lambda wc: expand(["results/filtered/{control}.sorted.bam"],
            control = get_controls_of_group(wc.group)),
    output:
        peak = "results/genrich/{group}.narrowPeak",
        bed = "results/genrich/{group}.bed",
    log:
        "logs/genrich/{group}.log"
    params:
        samples = lambda wc, input: ",".join(input.samples),
        controls = lambda wc, input: ",".join(input.controls),
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
        "Genrich {params.base} {params.single} {params.excl} {params.p_value} {params.q_value} -t {params.samples} -c {params.controls} -o {output.peak} 2>{log}"

 #The bedgraph output above has to be processed into 4 columns (first 3 + a datavalue of choice) for visualization  

#Only Frip score of treatment samples, ignoring controls or inputs (3 next rules).
rule bedtools_intersect:
    input:
        left="results/filtered/{sample}.sorted.bam",
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
        flagstats="results/filtered/stats/{sample}.mapped.flagstat"
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
        expand("results/bedtools_intersect/{sample}.intersected.bed",sample = samples.loc[samples['group'].isin(groups)].index)
    output:
        report("results/genrich/plots/plot_narrow_peaks_frip_score.pdf", caption="../report/plot_frip_score_genrich_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_narrow_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"

