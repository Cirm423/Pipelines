import random
import os

rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index=f"{assembly_path}star_genome_{assembly}",
    output:
        "results/star/pe/{samples_units}/Aligned.out.bam",
        "results/star/pe/{samples_units}/Aligned.toTranscriptome.out.bam",
        "results/star/pe/{samples_units}/SJ.out.tab",
        "results/star/pe/{samples_units}/Log.final.out",
    log:
        "logs/star-pe/{samples_units}.log",
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
        "results/star/se/{samples_units}/Aligned.out.bam",
        "results/star/se/{samples_units}/Aligned.toTranscriptome.out.bam",
        "results/star/se/{samples_units}/SJ.out.tab",
        "results/star/se/{samples_units}/Log.final.out",
    log:
        "logs/star-se/{samples_units}.log",
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
        sj=lambda wc: get_star_output_all_units(wc, fi='SJ',orig=True),
    output:
        "results/star/pe2/{samples_units}/Aligned.out.bam",
        "results/star/pe2/{samples_units}/Aligned.toTranscriptome.out.bam",
        "results/star/pe2/{samples_units}/Log.final.out",
    log:
        "logs/star-pe2/{samples_units}.log",
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
        sj=lambda wc: get_star_output_all_units(wc, fi='SJ',orig=True),
    output:
        "results/star/se2/{samples_units}/Aligned.out.bam",
        "results/star/se2/{samples_units}/Aligned.toTranscriptome.out.bam",
        "results/star/se2/{samples_units}/Log.final.out",
    log:
        "logs/star-se2/{samples_units}.log",
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input:"--quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} --sjdbFileChrStartEnd {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", input.sj, config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"

rule samtools_sort_pe:
    input:
        lambda wc: get_star_bam_uns(wc,original=True),
    output:
        "results/star/pe/{samples_units}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_sort_se:
    input:
        lambda wc: get_star_bam_uns(wc,original=True),
    output:
        "results/star/se/{samples_units}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
        tmp_dir = ""
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_sort_pe_2:
    input:
        get_star_bam_uns,
    output:
        "results/star/pe2/{samples_units}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
        tmp_dir = ""
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_sort_se_2:
    input:
        get_star_bam_uns,
    output:
        "results/star/se2/{samples_units}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
        tmp_dir = ""
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_index_pe:
    input:
        lambda wc: get_star_bam(wc, original=True)
    output:
        "results/star/pe/{samples_units}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/pe/{samples_units}.log",
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

rule samtools_index_se:
    input:
        lambda wc: get_star_bam(wc, original=True)
    output:
        "results/star/se/{samples_units}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/se/{samples_units}.log",
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

rule samtools_index_pe_2:
    input:
        get_star_bam
    output:
        "results/star/pe2/{samples_units}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/pe2/{samples_units}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

rule samtools_index_se_2:
    input:
        get_star_bam
    output:
        "results/star/se2/{samples_units}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/se2/{samples_units}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

rule rsem_pe:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        bam=get_star_transcript_bam,
        # bam = expand(
        #         "results/star/pe/{sample}-{unit}/Aligned.out.bam",
        #         unit=units["unit_name"],
        #         sample=units["sample_name"],
        #     ),
        # one of the index files created by rsem-prepare-reference; the file suffix is stripped and passed on to rsem
        #bai=get_star_bam_bai,
        reference=f"{assembly_path}rsem_reference/{assembly}.seq",
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results="results/rsem/pe¿/{samples_units}/mapped.genes.results",
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/rsem/pe¿/{samples_units}/mapped.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=True,
        out_path = lambda wildcards, output: os.path.dirname(output.genes_results) + "/" + "mapped",
        rsem_ref= lambda wildcards, input: os.path.splitext(input.reference)[0],
        # additional optional parameters to pass to rsem, for example,
        extra=f"--bam --seed {random.randint(0,100000)} --forward-prob {float(get_strandedness(units)[0])} {config['params']['rsem']}",
    threads: 24
    log:
        "logs/rsem/calculate_expression/{samples_units}.log",
    conda:
        "../envs/rsem.yaml"
    shell:
        "rsem-calculate-expression --num-threads {threads} {params.extra} --paired-end --alignments {input.bam} {params.rsem_ref} {params.out_path} > {log} 2>&1"

rule rsem_se:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        bam=get_star_transcript_bam,
        # bam= expand(
        #         "results/star/se/{sample}-{unit}/Aligned.out.bam",
        #         unit=units["unit_name"],
        #         sample=samples["sample_name"],
        #     ),
        # one of the index files created by rsem-prepare-reference; the file suffix is stripped and passed on to rsem
        #bai=get_star_bam_bai,
        reference=f"{assembly_path}rsem_reference/{assembly}.seq",
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results="results/rsem/se¿/{samples_units}/mapped.genes.results",
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/rsem/se¿/{samples_units}/mapped.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=False,
        out_path = lambda wildcards, output: os.path.dirname(output.genes_results) + "/" + "mapped",
        rsem_ref= lambda wildcards, input: os.path.splitext(input.reference)[0],
        # additional optional parameters to pass to rsem, for example,
        extra=f"--bam --seed {random.randint(0,100000)} --forward-prob {float(get_strandedness(units)[0])} {config['params']['rsem']}",
    threads: 24
    log:
        "logs/rsem/calculate_expression/{samples_units}.log",
    conda:
        "../envs/rsem.yaml"
    shell: 
        "rsem-calculate-expression --num-threads {threads} {params.extra} --alignments {input.bam} {params.rsem_ref} {params.out_path} > {log} 2>&1"