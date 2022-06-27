if genecode_assembly:

    rule get_genome_gencode:
        output:
            multiext(f"{assembly_path}{assembly}",".fa",".annotation.gtf"),
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
            f"{assembly_path}{assembly}.fa",
        output:
            f"{assembly_path}{assembly}.fa.fai",
        log:
            f"logs/genome-faidx_{config['resources']['ref']['assembly']}.log",
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
            f"logs/get-genome_{config['resources']['ref']['assembly']}.log",
        params:
            provider="UCSC",
            assembly=f"{config['resources']['ref']['assembly']}",
        cache: True
        conda:
            "../envs/genomepy.yaml"
        script:
            "../scripts/genomepy.py"


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
        f"logs/bwa_index_{config['resources']['ref']['assembly']}.log",
    params:
        algorithm="bwtsw",
    resources:
        mem_mb=369000,
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
        blacklist_path=f"{assembly_path}{assembly}.blacklist.bed"
    params:
        build = assembly,
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
        in_file=f"{assembly_path}{assembly}.blacklist.bed"
    output:
        f"{assembly_path}{assembly}.blacklist.sorted"
    params:
        extra=""
    log:
        "logs/ref/blacklist.sorted.log"
    wrapper:
        "v1.3.1/bio/bedtools/sort"

rule bedtools_format_blacklist:
    input:
        f"{assembly_path}{assembly}.blacklist.sorted"
    output:
        f"{assembly_path}{assembly}.blacklist_formated.sorted"
    log:
        "logs/ref/blacklist.format.log"
    run:
        transform=[]
        with open(input[0],"r") as f:
            for line in f:
                transform.append(line)
        
        with open(output[0],"w") as f:
            for line in transform:
                if line.startswith("chr"):
                    f.write(line)
                else:
                    if line.startswith("MT"):
                        line.replace("MT","chrM")
                    else:
                        f.write("chr" + line)

rule bedtools_complement_blacklist:
    input:
        in_file=f"{assembly_path}{assembly}.blacklist_formated.sorted",
        genome=f"{assembly_path}{assembly}.chrom.sizes.bedtools"
    output:
        f"{assembly_path}{assembly}.blacklist.sorted.complement"
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
        f"{assembly_path}{assembly}.gsize.txt"
    params:
        extra=config["resources"]["ref"]["macs-gsize"],
        build=assembly
    log:
        "logs/ref/gsize.log"
    conda:
        ""
    script:
        "../scripts/get_gsize.py"
