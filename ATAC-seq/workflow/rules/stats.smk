rule samtools_index:
    input:
        "results/mapped/{samples_units}.bam"
    output:
        "results/mapped/{samples_units}.bam.bai"
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/mapped/{samples_units}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"
        
rule samtools_flagstat:
    input:
        "results/mapped/{samples_units}.bam"
    output:
        "results/mapped/stats/{samples_units}.mapped.flagstat"
    log:
        "logs/samtools-flagstat/mapped/{samples_units}.mapped.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/mapped/{samples_units}.bam",
        idx = "results/mapped/{samples_units}.bam.bai"
    output:
        "results/mapped/stats/{samples_units}.mapped.idxstats"
    log:
        "logs/samtools-idxstats/mapped/{samples_units}.mapped.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/mapped/{samples_units}.bam"
    output:
        "results/mapped/stats/{samples_units}.mapped.stats.txt"
    params:
        extra=""
    log:
        "logs/samtools-stats/mapped/{samples_units}.mapped.log"
    wrapper:
        "v1.3.1/bio/samtools/stats"


rule samtools_index_merged:
    input:
        "results/merged_group/{group}.bam"
    output:
        "results/merged_group/{group}.bam.bai"
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/merged_group/{group}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"
        
rule samtools_flagstat_merged:
    input:
        "results/merged_group/{group}.bam"
    output:
        "results/merged/stats/{group}.merged_group.flagstat"
    log:
        "logs/samtools-flagstat/merged/{group}.merged_group.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"