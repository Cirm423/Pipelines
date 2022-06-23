if genecode_assembly:

    rule get_genome_gencode:
        output:
            multiext(f"{config['resources']['path']}{config['resources']['ref']['assembly']}",".fa",".annotation.gtf"),
        log:
            f"logs/get-genome_{config['resources']['ref']['assembly']}.log",
        params:
            assembly=f"{config['resources']['ref']['assembly']}",
        cache: True
        run:
            shell(f"wget -O {output[0]}.gz {genecode[config['resources']['ref']['assembly']]['assembly']} && gzip -d {output[0]}.gz")
            shell(f"wget -O {output[1]}.gz {genecode[config['resources']['ref']['assembly']]['gtf']} && gzip -d {output[1]}.gz")
    
    rule genome_faidx:
        input:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
        output:
            f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa.fai",
        log:
            f"logs/genome-faidx_{config['resources']['ref']['assembly']}.log",
        params:
            extra="",  # optional params string
        cache: True
        wrapper:
            "v1.3.1/bio/samtools/faidx"

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
            multiext(f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}", ".fa", ".fa.fai", ".fa.sizes",".annotation.gtf",".annotation.bed")
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


    rule move_annotation_ucsc:
        input:
            gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.annotation.gtf",
            bed=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.annotation.bed",
            sizes=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa.sizes",
            fai=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa.fai",
            fa=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa",
        output:
            multiext(f"{config['resources']['path']}{config['resources']['ref']['assembly']}",".annotation.gtf",".annotation.bed",".chrom.sizes",".fa.fai",".fa")
        params:
            folder=f"{config['resources']['path']}{config['resources']['ref']['assembly']}"
        cache: True
        log:
            f"logs/unzip_annotation_{config['resources']['ref']['assembly']}.log"
        shell:
            "mv {input.gtf} {output[0]} 2>{log} && mv {input.bed} {output[1]} 2>>{log} && mv {input.sizes} {output[2]} && mv {input.fai} {output[3]} && mv {input.fa} {output[4]} && rm -r {params.folder}"

rule bwa_index_meth:
    input:
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa",
    output:
        idx=multiext((f"{config['resources']['path']}{config['resources']['ref']['assembly']}.fa.bwameth"), ".c2t" ".c2t.amb", ".c2t.ann", ".c2t.bwt", ".c2t.pac", ".c2t.sa"),
    log:
        f"logs/bwa_meth_index_{config['resources']['ref']['assembly']}.log",
    resources:
        mem_mb=369000,
    cache: True
    conda:
        "../envs/bwa_meth.yaml"
    shell:
        "bwameth.py index {input}"


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
            "v1.3.1/bio/ucsc/faToTwoBit"

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
            "v1.3.1/bio/ucsc/twoBitInfo"

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
        "resources/sra-pe-reads/{accession}_1.fastq.gz",
        "resources/sra-pe-reads/{accession}_2.fastq.gz",
    params:
        extra="--skip-technical"
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "v1.3.1/bio/sra-tools/fasterq-dump"

rule sra_get_fastq_se:
    output:
        "resources/sra-se-reads/{accession}.fastq.gz"
    params:
        extra="--skip-technical"
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "v1.3.1/bio/sra-tools/fasterq-dump"