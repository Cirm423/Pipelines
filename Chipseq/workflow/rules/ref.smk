if genecode_assembly:

    rule get_genome_gencode:
        output:
            fasta=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
            gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf",
        log:
            f"logs/get-genome_{config['resources']['ref']['assembly']}.log",
        params:
            assembly=f"{config['resources']['ref']['assembly']}",
        cache: True
        run:
            shell(f"wget -O {output.fasta}.gz {genecode[config['resources']['ref']['assembly']]['assembly']} && gzip -d {output.fasta}.gz")
            shell(f"wget -O {output.gtf}.gz {genecode[config['resources']['ref']['assembly']]['gtf']} && gzip -d {output.gtf}.gz")
    
    rule genome_faidx:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        output:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa.fai",
        log:
            f"logs/genome-faidx_{config['resources']['ref']['assembly']}.log",
        cache: True
        wrapper:
            "0.77.0/bio/samtools/faidx"

    rule annot_gtf2bed:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf",
        output:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.bed",
        log:
            "logs/annot_gtf2bed.log",
        cache: True
        conda:
            "../envs/bedops.yaml"
        shell:
            r"""awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }}' {input} | gtf2bed - > {output}"""

else:

    rule get_genome_ucsc:
        output:
            multiext(f"{config['resources']['path']}{config['resources']['ref']['assembly']}", ".fa", ".fa.fai", ".fa.sizes"),
            temp(f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf.gz"),
            temp(f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.bed.gz"),
        log:
            f"logs/get-genome_{config['resources']['ref']['assembly']}.log",
        params:
            provider="UCSC",
            assembly=f"{config['resources']['ref']['assembly']}",
        cache: True
        conda:
            "../envs/genomepy.yaml"
        script:
            "../scripts/genomepy.py"


    rule unzip_annotation_ucsc:
        input:
            gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf.gz",
            bed=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.bed.gz",
            sizes=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa.sizes",
        output:
            gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.gtf",
            bed=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.annotation.bed",
            sizes=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes",
        cache: True
        log:
            f"logs/unzip_annotation_{config['resources']['ref']['assembly']}.log"
        shell:
            "gzip -dc {input.gtf} > {output.gtf} 2>{log} && gzip -dc {input.bed} > {output.bed} 2>>{log} && mv {input.sizes} {output.sizes}"

rule bwa_index:
    input:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
    output:
        multiext((f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"logs/bwa_index_{config['resources']['ref']['assembly']}.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"


if genecode_assembly:

    rule faToTwoBit_fa:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        output:
            temp(f"{config['resources']['path']}{config['resources']['ref']['assembly']}.2bit"),
        log:
            "logs/browser/fa_to_2bit.log"
        params:
            "" # optional params string
        wrapper:
            "0.78.0/bio/ucsc/faToTwoBit"

    rule twoBitInfo:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.2bit"
        output:
            temp(f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes.tmp")
        log:
            "logs/browser/chrom.sizes.log"
        params:
            "" # optional params string
        wrapper:
            "0.78.0/bio/ucsc/twoBitInfo"

    rule twoBitInfo_sort:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes.tmp"
        output:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes"
        cache: True
        shell:
            "sort -k2rn {input} > {output}"


rule twoBitInfo_sort_tobedtools:
    input:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes"
    output:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes.bedtools"
    cache: True
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


# SRA-download
rule sra_get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/sra-pe-reads/{accession}_1.fastq",
        "resources/sra-pe-reads/{accession}_2.fastq",
    params:
        extra=""
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump"

rule sra_get_fastq_se:
    output:
        "resources/sra-se-reads/{accession}.fastq"
    params:
        extra=""
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump"

rule generate_igenomes:
    output:
        f"{config['resources']['path']}igenomes.yaml"
    params:
        igenomes_release = config["resources"]["ref"]["igenomes_release"]
    log:
        "logs/ref/igenomes.log"
    conda:
        ""
    script:
        "../scripts/generate_igenomes.py"

rule generate_igenomes_blacklist:
    input:
        f"{config['resources']['path']}igenomes.yaml"
    output:
        blacklist_path=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.bed"
    params:
        build = config["resources"]["ref"]["assembly"],
        chromosome = config["resources"]["ref"]["chromosome"],
        blacklist = config["resources"]["ref"]["blacklist"]
    log:
        "logs/ref/blacklist.log"
    conda:
        ""
    script:
        "../scripts/generate_blacklist.py"

rule bedtools_sort_blacklist:
    input:
        in_file=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.bed"
    output:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.sorted"
    params:
        extra=""
    log:
        "logs/ref/blacklist.sorted.log"
    wrapper:
        "0.68.0/bio/bedtools/sort"

rule bedtools_complement_blacklist:
    input:
        in_file=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.sorted",
        genome=f"{config['resources']['path']}{config['resources']['ref']['assembly']}.chrom.sizes.bedtools"
    output:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.sorted.complement"
    params:
        extra=""
    log:
        "logs/ref/blacklist.sorted.complement.log"
    wrapper:
        "0.68.0/bio/bedtools/complement"

checkpoint get_gsize:
    input:
        f"{config['resources']['path']}igenomes.yaml"
    output:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.gsize.txt"
    params:
        extra=config["resources"]["ref"]["macs-gsize"],
        build=config["resources"]["ref"]["assembly"]
    log:
        "logs/ref/gsize.log"
    conda:
        ""
    script:
        "../scripts/get_gsize.py"
