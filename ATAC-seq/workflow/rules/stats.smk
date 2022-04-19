rule samtools_index:
    input:
        "results/filtered/{sample}.sorted.bam"
    output:
        "results/filtered/{sample}.sorted.bam.bai"
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/mapped/{sample}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"
        
rule samtools_flagstat:
    input:
        "results/filtered/{sample}.sorted.bam"
    output:
        "results/filtered/stats/{sample}.mapped.flagstat"
    log:
        "logs/samtools-flagstat/mapped/{sample}.mapped.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/filtered/{sample}.sorted.bam",
        idx = "results/filtered/{sample}.sorted.bam.bai"
    output:
        "results/filtered/stats/{sample}.mapped.idxstats"
    log:
        "logs/samtools-idxstats/mapped/{sample}.mapped.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/filtered/{sample}.sorted.bam"
    output:
        "results/filtered/stats/{sample}.mapped.stats.txt"
    params:
        extra=""
    log:
        "logs/samtools-stats/mapped/{sample}.mapped.log"
    wrapper:
        "v1.3.1/bio/samtools/stats"