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
            shell(f"wget -O {output[0]}.gz {genecode[config['ref']['assembly']]['assembly']} && gzip -d {output[0]}.gz")
            shell(f"wget -O {output[1]}.gz {genecode[config['ref']['assembly']]['gtf']} && gzip -d {output[1]}.gz")
    
    rule genome_faidx:
        input:
            f"{assembly_path}{assembly}.fa",
        output:
            f"{assembly_path}{assembly}.fa.fai",
        log:
            f"logs/genome-faidx_{assembly}.log",
        cache: True
        wrapper:
            "0.77.0/bio/samtools/faidx"

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
    #         gtf=f"{config['resources']}{assembly}/{assembly}.annotation.gtf",
    #         bed=f"{config['resources']}{assembly}/{assembly}.annotation.bed",
    #         sizes=f"{config['resources']}{assembly}/{assembly}.fa.sizes",
    #         fai=f"{config['resources']}{assembly}/{assembly}.fa.fai",
    #         fa=f"{config['resources']}{assembly}/{assembly}.fa",
    #     output:
    #         multiext(f"{config['resources']}{assembly}",".annotation.gtf",".annotation.bed",".chrom.sizes",".fa.fai",".fa")
    #     params:
    #         folder=f"{config['resources']}{assembly}"
    #     cache: True
    #     log:
    #         f"logs/unzip_annotation_{assembly}.log"
    #     shell:
    #         "mv {input.gtf} {output[0]} 2>{log} && mv {input.bed} {output[1]} 2>>{log} && mv {input.sizes} {output[2]} && mv {input.fai} {output[3]} && mv {input.fa} {output[4]} && rm -r {params.folder}"

rule bwa_index:
    input:
        f"{assembly_path}{assembly}.fa",
    output:
        multiext((f"{assembly_path}{assembly}.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"logs/bwa_index_{assembly}.log",
    threads: 12
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"


rule star_index:
    input:
        fasta=f"{assembly_path}{assembly}.fa",
        annotation=f"{assembly_path}{assembly}.annotation.gtf",
    output:
        directory(f"{assembly_path}star_genome_{assembly}"),
    threads: 24
    params:
        extra=f"--sjdbGTFfile {assembly_path}{assembly}.annotation.gtf --sjdbOverhang 100",
    log:
        f"logs/star_index_genome_{assembly}.log",
    cache: True
    wrapper:
        "0.77.0/bio/star/index"

rule rsem_ref:
    input:
        # reference FASTA with either the entire genome or transcript sequences
        f"{assembly_path}{assembly}.annotation.gtf",
        reference_genome=f"{assembly_path}{assembly}.fa",
    output:
        # one of the index files created and used by RSEM (required)
        multiext(f"{assembly_path}rsem_reference/{assembly}",".seq",".grp",".ti")
        # RSEM produces a number of other files which may optionally be specified as output (later 2 above); these may be provided so that snakemake is aware of them, but the wrapper doesn't do anything with this information other than to verify that the file path prefixes match that of output.seq.
    threads: 4
    params:
        # optional additional parameters, for example,
        #extra="--gtf annotations.gtf",
        # if building the index against a reference transcript set
        extra=f"--gtf {assembly_path}{assembly}.annotation.gtf",
        out_ref = f"{assembly_path}rsem_reference/{assembly}",
    log:
        f"logs/rsem/prepare-reference_{assembly}.log",
    conda:
        "../envs/rsem.yaml"
    cache: True
    shell:
        "rsem-prepare-reference --num-threads {threads} {params.extra} {input.reference_genome} {params.out_ref} > {log} 2>&1"

rule get_rmsk:
    output:
        temp(f"{assembly_path}{assembly}.rmsk.txt"),
    log:
        f"logs/get_rmsk_{assembly}.log",
    params:
        assembly = get_assembly_rmsk,
        dir = f"{assembly_path}",
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.assembly}/database/rmsk.txt.gz -P {params.dir} && gzip -dc {params.dir}rmsk.txt.gz > {output} && rm {params.dir}rmsk.txt.gz 2>{log}"

rule rmsk_to_bed:
    input:
        f"{assembly_path}{assembly}.rmsk.txt",
    output:
        f"{assembly_path}{assembly}.rmsk.bed",
    params:
        assemb={assembly}
    log:
        f"logs/rsem/rmsk_to_bed-{assembly}",
    cache: True
    shell:
        r"""awk 'BEGIN{{OFS="\t"}} {{print $6,$7,$8,$6";"$7";"$8";"$11";"$12";"$13,$2,$10}}' {input} > {output} 2>{log}"""

#Shell line for the rule bellow from here https://www.biostars.org/p/206342/

if genecode_assembly:

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


rule rmsk_bed2gtf:
    input:
        f"{assembly_path}{assembly}.rmsk.bed",
    output:
        f"{assembly_path}{assembly}.rmsk.gtf",
    log:
        "logs/rmsk_bed2gtf.log",
    cache: True
    conda:
        "../envs/ucscutils.yaml"
    shell:
        "bedToGenePred {input} stdout | genePredToGtf file stdin {output}"