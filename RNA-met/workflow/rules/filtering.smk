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

rule samtools_sort_genomic:
    input:
        "results/star/{star_lib}/{sample}/Aligned.out.bam"
    output:
        #temp("results/filtered/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
        "results/filtered/{star_lib}/{sample}.sortedByCoord.out.bam",
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"
rule remove_optical_duplicates_transcriptome:
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

rule remove_optical_duplicates:
    input:
        bams="results/filtered/{star_lib}/{sample}.sortedByCoord.out.bam"
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
