rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        idx=f"{assembly_path}star_genome_{assembly}",
    output:
        aln=temp("results/star/pe/{sample}/Aligned.out.bam"),
        aln_trans=temp("results/star/pe/{sample}/Aligned.toTranscriptome.out.bam"),
        sj="results/star/pe/{sample}/SJ.out.tab",
        log="results/star/pe/{sample}/Log.out",
        log_final="results/star/pe/{sample}/Log.final.out",
    log:
        "logs/star-pe/{sample}.log",
    params:
        extra="--outSAMunmapped Within KeepPairs --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", config["params"]["star"]
        ),
    threads: 24
    resources:
        mem_mb = 5000
    conda:
        "../envs/star.yaml"
    script:
        "../scripts/start.py"

rule align_se:
    input:
        fq1=get_map_reads_input_R1,
        idx=f"{assembly_path}star_genome_{assembly}",
    output:
        aln=temp("results/star/se/{sample}/Aligned.out.bam"),
        aln_trans=temp("results/star/se/{sample}/Aligned.toTranscriptome.out.bam"),
        sj="results/star/se/{sample}/SJ.out.tab",
        log="results/star/se/{sample}/Log.out",
        log_final="results/star/se/{sample}/Log.final.out",
    log:
        "logs/star-se/{sample}.log",
    params:
        extra="--quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --sjdbGTFfile {} {}".format(
            f"{assembly_path}{assembly}.annotation.gtf", config["params"]["star"]
        ),
    threads: 24
    resources:
        mem_mb = 5000
    conda:
        "../envs/star.yaml"
    script:
        "../scripts/start.py"

rule samtools_sort_star:
    input:
        "results/star/{lib}/{sample}/Aligned.out.bam"
    output:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "v2.6.0/bio/samtools/sort"

rule samtools_index_star:
    input:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/{lib}/{sample}.log",
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"

