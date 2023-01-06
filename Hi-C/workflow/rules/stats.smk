rule samtools_flagstat:
    input:
        "results/mapped/{sample}_R{read}.sorted.bam"
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
        bam = "results/mapped/{sample}_R{read}.sorted.bam",
        idx = "results/mapped/{sample}_R{read}.sorted.bam.bai"
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
        "results/mapped/{sample}_R{read}.sorted.bam"
    output:
        "results/mapped/{sample}_R{read}.stats.txt"
    params:
        extra=""
    log:
        "logs/samtools-stats/{sample}_R{read}.log"
    wrapper:
        "v1.3.1/bio/samtools/stats"