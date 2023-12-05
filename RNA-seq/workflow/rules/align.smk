import random
import os

rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index=f"{assembly_path}star_genome_{assembly}",
    output:
        "results/star/pe/{sample}/Aligned.out.bam",
        "results/star/pe/{sample}/Aligned.toTranscriptome.out.bam",
        "results/star/pe/{sample}/SJ.out.tab",
        "results/star/pe/{sample}/Log.final.out",
    log:
        "logs/star-pe/{sample}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--outSAMunmapped Within KeepPairs --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"


rule align_se:
    input:
        fq1=get_map_reads_input_R1,
        index=f"{assembly_path}star_genome_{assembly}",
    output:
        "results/star/se/{sample}/Aligned.out.bam",
        "results/star/se/{sample}/Aligned.toTranscriptome.out.bam",
        "results/star/se/{sample}/SJ.out.tab",
        "results/star/se/{sample}/Log.final.out",
    log:
        "logs/star-se/{sample}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"

rule align_pe_2pass:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index=f"{assembly_path}star_genome_{assembly}",
        sj=expand("results/star/pe/{sample}/SJ.out.tab",sample=samples.sample_name)
    output:
        "results/star/pe2/{sample}/Aligned.out.bam",
        "results/star/pe2/{sample}/Aligned.toTranscriptome.out.bam",
        "results/star/pe2/{sample}/Log.final.out",
    log:
        "logs/star-pe2/{sample}.log",
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input:"--quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} --sjdbFileChrStartEnd {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", input.sj, config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"


rule align_se_2pass:
    input:
        fq1=get_map_reads_input_R1,
        index=f"{assembly_path}star_genome_{assembly}",
        sj=expand("results/star/se/{sample}/SJ.out.tab",sample=samples.sample_name)
    output:
        "results/star/se2/{sample}/Aligned.out.bam",
        "results/star/se2/{sample}/Aligned.toTranscriptome.out.bam",
        "results/star/se2/{sample}/Log.final.out",
    log:
        "logs/star-se2/{sample}.log",
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input:"--quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} --sjdbFileChrStartEnd {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", input.sj, config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"

rule samtools_sort_star:
    input:
        "results/star/{lib}/{sample}/Aligned.out.bam"
    output:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_index_star:
    input:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/{lib}/{sample}.log",
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

rule rsem:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        bam="results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.fixed.out.bam",
        # bam = expand(
        #         "results/star/pe/{sample}-{unit}/Aligned.out.bam",
        #         unit=units["unit_name"],
        #         sample=units["sample_name"],
        #     ),
        # one of the index files created by rsem-prepare-reference; the file suffix is stripped and passed on to rsem
        reference=f"{assembly_path}rsem_reference/{assembly}.seq",
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results="results/rsem/{star_lib}/{sample}/mapped.genes.results",
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/rsem/{star_lib}/{sample}/mapped.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired_end = not config["single_end"],
        out_path = lambda wildcards, output: os.path.dirname(output.genes_results) + "/" + "mapped",
        rsem_ref= lambda wildcards, input: os.path.splitext(input.reference)[0],
        # additional optional parameters to pass to rsem, for example,
        extra=f"--bam --forward-prob {float(get_strandedness(units)[0])} {config['params']['rsem']}",
    threads: 24
    log:
        "logs/rsem/calculate_expression/{sample}-{star_lib}.log",
    conda:
        "../envs/rsem.yaml"
    shell:
        "rsem-calculate-expression --num-threads {threads} {params.extra} --paired-end --alignments {input.bam} {params.rsem_ref} {params.out_path} > {log} 2>&1"
