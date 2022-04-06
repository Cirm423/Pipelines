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