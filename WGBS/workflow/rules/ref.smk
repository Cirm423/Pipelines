if genecode_assembly:

    rule get_genome_gencode:
        output:
            multiext(f"{assembly_path}{assembly}", ".fa", ".annotation.gtf"),
        log:
            f"logs/get-genome_{assembly}.log",
        params:
            assembly=f"{assembly}",
        cache: True
        run:
            shell(
                f"wget -O {output[0]}.gz {genecode[assembly]['assembly']} && gzip -d {output[0]}.gz"
            )
            shell(
                f"wget -O {output[1]}.gz {genecode[assembly]['gtf']} && gzip -d {output[1]}.gz"
            )
            shell(f"chmod -R g+w {assembly_path} || true")

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
            multiext(
                f"{assembly_path}{assembly}",
                ".fa",
                ".fa.fai",
                ".fa.sizes",
                ".annotation.gtf",
                ".annotation.bed",
            ),
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
            f"{assembly_path}{assembly}.fa.sizes",
        output:
            f"{assembly_path}{assembly}.chrom.sizes",
        cache: True
        shell:
            "mv {input} {output}"


# rule move_annotation_ucsc:
#     input:
#         gtf=f"{assembly_path}{assembly}.annotation.gtf",
#         bed=f"{assembly_path}{assembly}.annotation.bed",
#         sizes=f"{assembly_path}{assembly}.fa.sizes",
#         fai=f"{assembly_path}{assembly}.fa.fai",
#         fa=f"{assembly_path}{assembly}.fa",
#     output:
#         multiext(f"{config['resources']['path']}{assembly}",".annotation.gtf",".annotation.bed",".chrom.sizes",".fa.fai",".fa")
#     params:
#         folder=f"{config['resources']['path']}{assembly}"
#     cache: True
#     log:
#         f"logs/unzip_annotation_{assembly}.log"
#     shell:
#         "mv {input.gtf} {output[0]} 2>{log} && mv {input.bed} {output[1]} 2>>{log} && mv {input.sizes} {output[2]} && mv {input.fai} {output[3]} && mv {input.fa} {output[4]} && rm -r {params.folder}"


rule get_phage_genome:
    output:
        f"{phage_path}phage_lambda.fa",
    cache: True
    log:
        "get_lamda_genome.log",
    shell:
        "wget https://www.ebi.ac.uk/ena/browser/api/fasta/J02459.1 -O {output}"


rule bwa_index_meth:
    input:
        f"{assembly_path}{assembly}.fa",
    output:
        idx=multiext(
            (f"{assembly_path}{assembly}.fa.bwameth"),
            ".c2t",
            ".c2t.amb",
            ".c2t.ann",
            ".c2t.bwt",
            ".c2t.pac",
            ".c2t.sa",
        ),
    log:
        f"logs/bwa_meth_index_{assembly}.log",
    threads: 24
    cache: True
    conda:
        "../envs/bwa_meth.yaml"
    shell:
        "bwameth.py index {input} 2>{log}"


rule bismark_genome_preparation_fa:
    input:
        f"{assembly_path}{assembly}.fa",
    output:
        directory(f"{assembly_path}Bisulfite_Genome"),
    log:
        f"logs/bismark/{assembly}_index.log",
    params:
        "",  # optional params string
    threads: 24
    cache: True
    wrapper:
        "v2.2.0/bio/bismark/bismark_genome_preparation"


rule bwa_index_meth_phage:
    input:
        f"{phage_path}phage_lambda.fa",
    output:
        idx=multiext(
            (f"{phage_path}phage_lambda.fa.bwameth"),
            ".c2t",
            ".c2t.amb",
            ".c2t.ann",
            ".c2t.bwt",
            ".c2t.pac",
            ".c2t.sa",
        ),
    log:
        f"logs/bwa_meth_index_phage.log",
    threads: 1
    cache: True
    conda:
        "../envs/bwa_meth.yaml"
    shell:
        "bwameth.py index {input} 2>{log}"


rule bismark_genome_preparation_phage:
    input:
        f"{phage_path}phage_lambda.fa",
    output:
        directory(f"{phage_path}Bisulfite_Genome"),
    log:
        f"logs/bismark/phage_index.log",
    params:
        "",  # optional params string
    threads: 1
    cache: True
    wrapper:
        "v2.2.0/bio/bismark/bismark_genome_preparation"


if genecode_assembly:

    rule faToTwoBit_fa:
        input:
            f"{assembly_path}{assembly}.fa",
        output:
            temp(f"{assembly_path}{assembly}.2bit"),
        log:
            "logs/browser/fa_to_2bit.log",
        params:
            "",  # optional params string
        wrapper:
            "v1.3.1/bio/ucsc/faToTwoBit"

    rule twoBitInfo:
        input:
            f"{assembly_path}{assembly}.2bit",
        output:
            temp(f"{assembly_path}{assembly}.chrom.sizes.tmp"),
        log:
            "logs/browser/chrom.sizes.log",
        params:
            "",  # optional params string
        wrapper:
            "v1.3.1/bio/ucsc/twoBitInfo"

    rule twoBitInfo_sort:
        input:
            f"{assembly_path}{assembly}.chrom.sizes.tmp",
        output:
            f"{assembly_path}{assembly}.chrom.sizes",
        cache: True
        shell:
            "sort -k2rn {input} > {output}"


rule twoBitInfo_sort_tobedtools:
    input:
        f"{assembly_path}{assembly}.chrom.sizes",
    output:
        f"{assembly_path}{assembly}.chrom.sizes.bedtools",
    cache: True
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


# Make the bed12 necessary for genomation
rule gtf2bed12:
    input:
        f"{assembly_path}{assembly}.annotation.gtf",
    output:
        f"{assembly_path}{assembly}.annotation.bed12",
    log:
        "logs/gtf2bed12.log",
    cache: True
    conda:
        "../envs/ucscutils.yaml"
    shell:
        "gtfToGenePred {input} stdout | genePredToBed stdin {output}"


# SRA-download
rule sra_get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/sra-pe-reads/{accession}_1.fastq.gz",
        "resources/sra-pe-reads/{accession}_2.fastq.gz",
    params:
        extra="--skip-technical",
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log",
    wrapper:
        "v1.3.1/bio/sra-tools/fasterq-dump"


rule sra_get_fastq_se:
    output:
        "resources/sra-se-reads/{accession}.fastq.gz",
    params:
        extra="--skip-technical",
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log",
    wrapper:
        "v1.3.1/bio/sra-tools/fasterq-dump"
