if genecode_assembly:

    rule get_genome_gencode:
        output:
            fasta=f"{config['resources']}{config['ref']['assembly']}.fa",
            gtf=f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
        log:
            f"logs/get-genome_{config['ref']['assembly']}.log",
        params:
            assembly=f"{config['ref']['assembly']}",
        cache: True
        run:
            shell(f"wget -O {output.fasta}.gz {genecode[config['ref']['assembly']]['assembly']} && gzip -d {output.fasta}.gz")
            shell(f"wget -O {output.gtf}.gz {genecode[config['ref']['assembly']]['gtf']} && gzip -d {output.gtf}.gz")
    
    rule genome_faidx:
        input:
            f"{config['resources']}{config['ref']['assembly']}.fa",
        output:
            f"{config['resources']}{config['ref']['assembly']}.fa.fai",
        log:
            f"logs/genome-faidx_{config['ref']['assembly']}.log",
        cache: True
        wrapper:
            "0.77.0/bio/samtools/faidx"

else:

    rule get_genome_ucsc:
        output:
            multiext(f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}", ".fa", ".fa.fai", ".fa.sizes"),
            temp(f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.annotation.gtf.gz"),
            temp(f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.annotation.bed.gz"),
        log:
            f"logs/get-genome_{config['ref']['assembly']}.log",
        params:
            provider="UCSC",
            assembly=f"{config['ref']['assembly']}",
        cache: True
        conda:
            "../envs/genomepy.yaml"
        script:
            "../scripts/genomepy.py"


    rule unzip_annotation_ucsc:
        input:
            gtf=f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.annotation.gtf.gz",
            bed=f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.annotation.bed.gz",
            sizes=f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.fa.sizes",
            fai=f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.fa.fai",
            fa=f"{config['resources']}{config['ref']['assembly']}/{config['ref']['assembly']}.fa",
        output:
            gtf=f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
            bed=f"{config['resources']}{config['ref']['assembly']}.annotation.bed",
            sizes=f"{config['resources']}{config['ref']['assembly']}.chrom.sizes",
            fai=f"{config['resources']}{config['ref']['assembly']}.fa.fai",
            fa=f"{config['resources']}{config['ref']['assembly']}.fa",
        cache: True
        params:
            folder=f"{config['resources']}{config['ref']['assembly']}"
        log:
            f"logs/unzip_annotation_{config['ref']['assembly']}.log"
        shell:
            "gzip -dc {input.gtf} > {output.gtf} 2>{log} && gzip -dc {input.bed} > {output.bed} 2>>{log} && mv {input.sizes} {output.sizes} && mv {input.fai} {output.fai} && mv {input.fa} {output.fa} && rm -r {params.folder}"

rule bwa_index:
    input:
        f"{config['resources']}{config['ref']['assembly']}.fa",
    output:
        multiext((f"{config['resources']}{config['ref']['assembly']}.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"logs/bwa_index_{config['ref']['assembly']}.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"


rule star_index:
    input:
        fasta=f"{config['resources']}{config['ref']['assembly']}.fa",
        annotation=f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
    output:
        directory(f"{config['resources']}star_genome_{config['ref']['assembly']}"),
    threads: 24
    params:
        extra=f"--sjdbGTFfile {config['resources']}{config['ref']['assembly']}.annotation.gtf --sjdbOverhang 100",
    log:
        f"logs/star_index_genome_{config['ref']['assembly']}.log",
    cache: True
    wrapper:
        "0.77.0/bio/star/index"

rule rsem_ref:
    input:
        # reference FASTA with either the entire genome or transcript sequences
        f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
        reference_genome=f"{config['resources']}{config['ref']['assembly']}.fa",
    output:
        # one of the index files created and used by RSEM (required)
        seq=f"{config['resources']}rsem_reference_{config['ref']['assembly']}.seq",
        # RSEM produces a number of other files which may optionally be specified as output; these may be provided so that snakemake is aware of them, but the wrapper doesn't do anything with this information other than to verify that the file path prefixes match that of output.seq.
        # for example,
        grp=f"{config['resources']}rsem_reference_{config['ref']['assembly']}.grp",
        ti=f"{config['resources']}rsem_reference_{config['ref']['assembly']}.ti",
    threads: 4
    params:
        # optional additional parameters, for example,
        #extra="--gtf annotations.gtf",
        # if building the index against a reference transcript set
        extra=f"--gtf {config['resources']}{config['ref']['assembly']}.annotation.gtf",
        out_ref = f"{config['resources']}rsem_reference_{config['ref']['assembly']}",
    log:
        f"logs/rsem/prepare-reference_{config['ref']['assembly']}.log",
    cache: True
    shell:
        "rsem-prepare-reference --num-threads {threads} {params.extra} {input.reference_genome} {params.out_ref} > {log} 2>&1"

rule get_rmsk:
    output:
        temp(f"{config['resources']}{config['ref']['assembly']}.rmsk.txt"),
    log:
        f"logs/get_rmsk_{config['ref']['assembly']}.log",
    params:
        assembly = get_assembly_rmsk,
        dir = f"{config['resources']}",
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.assembly}/database/rmsk.txt.gz -P {params.dir} && gzip -dc {params.dir}rmsk.txt.gz > {output} && rm {params.dir}rmsk.txt.gz 2>{log}"

rule rmsk_to_bed:
    input:
        f"{config['resources']}{config['ref']['assembly']}.rmsk.txt",
    output:
        f"{config['resources']}{config['ref']['assembly']}.rmsk.bed",
    params:
        assemb={config['ref']['assembly']}
    log:
        f"logs/rsem/rmsk_to_bed-{config['ref']['assembly']}",
    cache: True
    shell:
        r"""awk 'BEGIN{{OFS="\t"}} {{print $6,$7,$8,$6";"$7";"$8";"$11";"$12";"$13,$2,$10}}' {input} > {output} 2>{log}"""

#Shell line for the rule bellow from here https://www.biostars.org/p/206342/

if genecode_assembly:

    rule annot_gtf2bed:
        input:
            f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
        output:
            f"{config['resources']}{config['ref']['assembly']}.annotation.bed",
        log:
            "logs/annot_gtf2bed.log",
        cache: True
        conda:
            "../envs/bedops.yaml"
        shell:
            r"""awk '{{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }}' {input} | gtf2bed - > {output}"""


rule rmsk_bed2gtf:
    input:
        f"{config['resources']}{config['ref']['assembly']}.rmsk.bed",
    output:
        f"{config['resources']}{config['ref']['assembly']}.rmsk.gtf",
    log:
        "logs/rmsk_bed2gtf.log",
    cache: True
    conda:
        "../envs/ucscutils.yaml"
    shell:
        "bedToGenePred {input} stdout | genePredToGtf file stdin {output}"