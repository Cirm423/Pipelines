rule bwa_mem:
    input:
        reads = get_map_reads_input if config["params"]["trimming"]["activate"] else get_fastqs,
        idx = rules.bwa_index.output
    output:
        temp("results/mapped/{sample}-{unit}.bam")
    log:
        "logs/bwa/bwa_mem/{sample}-{unit}.log"
    params:
        #index= lambda w, input: os.path.splitext(input.idx[0])[0],
        extra= get_read_group,
        sorting="samtools",
        sort_order="coordinate",
        sort_extra="",
    threads: 8
    wrapper:
        "v1.21.1/bio/bwa/mem"

rule merge_bams:
    input:
        lambda w: expand("results/mapped/{sample}-{unit}.bam",
            sample = w.sample,
            unit = units.loc[units['sample'] == w.sample].unit.to_list()
        )
    output:
        temp("results/merged/{sample}.bam")
    log:
        "logs/picard/mergebamfiles/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate",
    wrapper:
        "v0.87.0/bio/picard/mergesamfiles"

rule mark_merged_duplicates:
    input:
        bams="results/merged/{sample}.bam"
    output:
        bam=temp("results/picard_dedup/{sample}.bam"),
        metrics="results/picard_dedup/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        extra=f"--REMOVE_DUPLICATES false --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE {config['params']['optical_distance']}",
    threads: 4
    wrapper:
        "v3.0.2/bio/picard/markduplicates"
