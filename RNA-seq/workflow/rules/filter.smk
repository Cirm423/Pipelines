rule samtools_sort_transcriptome:
    input:
        "results/star/{star_lib}/{sample}/Aligned.toTranscriptome.out.bam"
    output:
        #temp("results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
        "results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam",
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule remove_optical_duplicates:
    input:
        bams="results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"
    output:
        #bam=temp("results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
        bam=temp("results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByCoord.out.bam"),
        metrics=report("results/filtered/{star_lib}/{sample}.toTranscriptome.metrics.txt",category="Optical replicates")
    log:
        "logs/picard/filtered/{sample}_{star_lib}.log"
    params:
        extra=f"--REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE {config['params']['optical_distance']} --REMOVE_SEQUENCING_DUPLICATES true",
    threads: 4
    resources:
        mem_mb = 5000
    wrapper:
        "v3.0.2/bio/picard/markduplicates"

# rule samtools_view_filter:
#     input:
#         "results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"
#     output:
#         "results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByCoord.out.bam"
#     params:
#         extra="-e '![DT]'"
#     log:
#         "logs/samtools-view/{sample}_{star_lib}.log"
#     wrapper:
#         "v1.3.1/bio/samtools/view"

rule convert_bam_for_rsem:
    input:
        "results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByCoord.out.bam"
    output:
        temp("results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.fixed.out.bam"),
    params:
        out_file = lambda wildcards, output: os.path.splitext(output[0])[0],
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    conda:
        "../envs/rsem.yaml"
    shell:
        "convert-sam-for-rsem -p {threads} {input} {params.out_file}"


rule remove_optical_duplicates_TE:
    input:
        bams="results/star/{star_lib}/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        #bam=temp("results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
        bam=temp("results/filtered/{star_lib}/{sample}.filtered.sortedByCoord.out.bam"),
        metrics=report("results/filtered/{star_lib}/{sample}.metrics.txt",category="Optical replicates")
    log:
        "logs/picard/filtered/{sample}_{star_lib}.log"
    params:
        extra=f"--REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE {config['params']['optical_distance']} --REMOVE_SEQUENCING_DUPLICATES true",
    threads: 4
    resources:
        mem_mb = 5000
    wrapper:
        "v3.0.2/bio/picard/markduplicates"