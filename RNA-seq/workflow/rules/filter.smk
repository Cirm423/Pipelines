rule samtools_sort_transcriptome:
    input:
        "results/star/{star_lib}/{sample}/Aligned.toTranscriptome.out.bam"
    output:
        temp("results/picard_dedup/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule mark_optical_duplicates:
    input:
        "results/picard_dedup/{star_lib}/{sample}/Aligned.toTranscriptome.sortedByCoord.out.bam"
    output:
        bam=temp("results/picard_dedup/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"),
        metrics=report("results/picard_dedup/{star_lib}/{sample}.metrics.txt",category="Optical replicates")
    log:
        "logs/picard/picard_dedup/{sample}_{star_lib}.log"
    params:
        f"REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE={config['params']['optical_distance']} --TAGGING_POLICY=OpticalOnly",
    threads: 4
    wrapper:
        "v0.87.0/bio/picard/markduplicates"

rule samtools_view_filter:
    input:
        "results/picard_dedup/{star_lib}/{sample}.toTranscriptome.sortedByCoord.out.bam"
    output:
        "results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByCoord.out.bam"
    params:
        extra="-e '![DT]'"
    log:
        "logs/samtools-view/{sample}_{star_lib}.log"
    wrapper:
        "v1.3.1/bio/samtools/view"

rule samtools_sort_transcriptome_name:
    input:
        "results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByCoord.out.bam"
    output:
        "results/filtered/{star_lib}/{sample}.toTranscriptome.filtered.sortedByName.out.bam",
    params:
        extra = "-n",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"