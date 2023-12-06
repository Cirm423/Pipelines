import os

if genecode_assembly:

    rule faToTwoBit_fa:
        input:
            f"{assembly_path}{assembly}.fa",
        output:
            temp(f"{assembly_path}{assembly}.2bit"),
        log:
            "logs/browser/fa_to_2bit.log"
        params:
            "" # optional params string
        wrapper:
            "v2.6.0/bio/ucsc/faToTwoBit"

    rule twoBitInfo:
        input:
            f"{assembly_path}{assembly}.2bit"
        output:
            temp(f"{assembly_path}{assembly}.chrom.sizes.tmp")
        log:
            "logs/browser/chrom.sizes.log"
        params:
            "" # optional params string
        wrapper:
            "v2.6.0/bio/ucsc/twoBitInfo"

    rule twoBitInfo_sort:
        input:
            f"{assembly_path}{assembly}.chrom.sizes.tmp"
        output:
            f"{assembly_path}{assembly}.chrom.sizes"
        cache: True
        shell:
            "sort -k2rn {input} > {output}"

rule BamCoverage_str1:
    input:
        expand(["results/star/{star_lib}/{{sample}}/Aligned.sortedByCoord.out.bam","results/star/{star_lib}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"],
            star_lib=star_lib)
    output:
        "results/browser/{sample}.str1.bw",
    params:
        norm = config["params"]["bamcoverage"],
        stranded = "--exactScaling" if get_strandedness(units)[0] == 0.5 else "--filterRNAstrand forward --exactScaling"
    log: 
        "logs/browser/{sample}.BamCoverage.log"
    conda:
        "../envs/deeptools.yaml"
    threads: 12
    shell:
        "bamCoverage -b {input[0]} -o {output} -of bigwig -p {threads} {params.norm} {params.stranded} 2>{log}"

rule BamCoverage_str2:
    input:
        expand(["results/star/{star_lib}/{{sample}}/Aligned.sortedByCoord.out.bam","results/star/{star_lib}/{{sample}}/Aligned.sortedByCoord.out.bam.bai"],
            star_lib=star_lib)
    output:
        "results/browser/{sample}.str2.bw",
    params:
        norm = config["params"]["bamcoverage"],
        stranded = "--filterRNAstrand reverse --exactScaling"
    log: 
        "logs/browser/{sample}.BamCoverage.log"
    conda:
        "../envs/deeptools.yaml"
    threads: 12
    shell:
        "bamCoverage -b {input[0]} -o {output} -of bigwig -p {threads} {params.norm} {params.stranded} 2>{log}"
