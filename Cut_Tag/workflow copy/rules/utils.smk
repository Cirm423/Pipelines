rule samtools_index:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step}/{samples_units}.bam.bai"
    params:
        extra="" # optional params string
    log:
        "logs/samtools-index/{step}/{samples_units}.log"
    wrapper:
        "v1.3.1/bio/samtools/index"
