if genecode_assembly:

    rule get_genome_gencode:
        output:
            multiext(f"{assembly_path}{assembly}",".fa",".annotation.gtf"),
        log:
            f"logs/get-genome_{assembly}.log",
        params:
            assembly=f"{assembly}",
        cache: True
        run:
            shell(f"wget -O {output[0]}.gz {genecode[assembly]['assembly']} && gzip -d {output[0]}.gz")
            shell(f"wget -O {output[1]}.gz {genecode[assembly]['gtf']} && gzip -d {output[1]}.gz")
            shell(f"chmod -R g+w {assembly_path}")
    
    rule genome_faidx:
        input:
            f"{assembly_path}{assembly}.fa",
        output:
            f"{assembly_path}{assembly}.fa.fai",
        log:
            f"logs/genome-faidx_{assembly}.log",
        params:
            extra="",  # optional params string
        cache: True
        wrapper:
            "v1.3.1/bio/samtools/faidx"

    rule annot_gtf2bed:
        input:
            f"{assembly_path}{assembly}.annotation.gtf",
        output:
            f"{assembly_path}{assembly}.annotation.bed",
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
            multiext(f"{assembly_path}{assembly}", ".fa", ".fa.fai", ".fa.sizes",".annotation.gtf",".annotation.bed")
        log:
            f"logs/get-genome_{assembly}.log",
        params:
            provider="UCSC",
            assembly=f"{assembly}",
        cache: True
        conda:
            "../envs/genomepy.yaml"
        script:
            "../scripts/genomepy.py"

    rule rename_sizes:
        input:
            f"{assembly_path}{assembly}.fa.sizes"
        output:
            f"{assembly_path}{assembly}.chrom.sizes"
        cache: True
        shell:
            "mv {input} {output}"

    # rule move_annotation_ucsc:
    #     input:
    #         gtf=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.annotation.gtf",
    #         bed=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.annotation.bed",
    #         sizes=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa.sizes",
    #         fai=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa.fai",
    #         fa=f"{config['resources']['path']}{config['resources']['ref']['assembly']}/{config['resources']['ref']['assembly']}.fa",
    #     output:
    #         multiext(f"{config['resources']['path']}{config['resources']['ref']['assembly']}",".annotation.gtf",".annotation.bed",".chrom.sizes",".fa.fai",".fa")
    #     params:
    #         folder=f"{config['resources']['path']}{config['resources']['ref']['assembly']}"
    #     cache: True
    #     log:
    #         f"logs/unzip_annotation_{config['resources']['ref']['assembly']}.log"
    #     shell:
    #         "mv {input.gtf} {output[0]} 2>{log} && mv {input.bed} {output[1]} 2>>{log} && mv {input.sizes} {output[2]} && mv {input.fai} {output[3]} && mv {input.fa} {output[4]} && rm -r {params.folder}"

rule bwa_index:
    input:
        f"{assembly_path}{assembly}.fa",
    output:
        idx=multiext((f"{assembly_path}{assembly}.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"logs/bwa_index_{assembly}.log",
    params:
        algorithm="bwtsw",
    threads: 12
    cache: True
    wrapper:
        "v1.3.1/bio/bwa/index"


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
            "v1.3.1/bio/ucsc/faToTwoBit"

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
            "v1.3.1/bio/ucsc/twoBitInfo"

    rule twoBitInfo_sort:
        input:
            f"{assembly_path}{assembly}.chrom.sizes.tmp"
        output:
            f"{assembly_path}{assembly}.chrom.sizes"
        cache: True
        shell:
            "sort -k2rn {input} > {output}"


rule twoBitInfo_sort_tobedtools:
    input:
        f"{assembly_path}{assembly}.chrom.sizes"
    output:
        f"{assembly_path}{assembly}.chrom.sizes.bedtools"
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

# Juicer jar file

rule get_juicer_jar:
    output:
        "resources/juicer/juicer_tools.2.20.00.jar"
    log:
        "logs/get_juicer_jar.log"
    shell:
        "wget https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/juicer_tools.2.20.00.jar -P ./resources/juicer"

#Not used in favor of fanc directionality implementation

#Install domaincaller in the env with the dependencies (only available in pip)

# rule domaincaller_init:
#     output:
#         temp(touch("results/domaincaller/package.done"))
#     log:
#         "logs/domaincaller/domaincaller_init.log"
#     conda:
#         "../envs/domaincaller.yaml"
#     shell:
#         "pip install domaincaller 2>{log}"