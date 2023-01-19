rule sort_mapped_stats:
    input:
        "results/mapped/{sample}_R{read}.bam"
    output:
        temp("results/mapped/{sample}_R{read}.sorted_stats.bam")
    params:
        extra=""
    log:
        "logs/bamtools_filtered/{sample}_R{read}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule samtools_flagstat:
    input:
        "results/mapped/{sample}_R{read}.sorted_stats.bam"
    output:
        "results/mapped/{sample}_R{read}.flagstat"
    log:
        "logs/samtools-flagstat/{sample}_R{read}.mapped.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/mapped/{sample}_R{read}.sorted_stats.bam",
        idx = "results/mapped/{sample}_R{read}.sorted_stats.bam.bai"
    output:
        "results/mapped/{sample}_R{read}.idxstats"
    log:
        "logs/samtools-idxstats/{sample}_R{read}.log"
    params:
        extra=""
    wrapper:
        "v1.3.1/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/mapped/{sample}_R{read}.sorted_stats.bam"
    output:
        "results/mapped/{sample}_R{read}.stats.txt"
    params:
        extra=""
    log:
        "logs/samtools-stats/{sample}_R{read}.log"
    wrapper:
        "v1.3.1/bio/samtools/stats"