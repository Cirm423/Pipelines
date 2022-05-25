rule samtools_flagstat:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.flagstat"
    log:
        "logs/samtools-flagstat/{step}/{samples_units}.{step}.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/{step}/{samples_units}.bam",
        idx = "results/{step}/{samples_units}.bam.bai"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.idxstats"
    log:
        "logs/samtools-idxstats/{step}/{samples_units}.{step}.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.stats.txt"
    params:
        extra=""
    log:
        "logs/samtools-stats/{step}/{samples_units}.{step}.log"
    wrapper:
        "v1.3.1/bio/samtools/stats"

rule samtools_flagstat_spike:
    input:
        "results/mapped/{sample}_spike-in.bam"
    output:
        "results/mapped/{sample}_spike-in.bam.flagstat"
    log:
        "logs/samtools-flagstat/mapped/{sample}_spike-in.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"