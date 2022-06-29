rule get_genome:
    output:
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta",
    log:
        f"logs/get-genome_{config['ref']['build']}_{config['ref']['release']}.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        f"logs/get_annotation_{config['ref']['build']}_{config['ref']['release']}.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta",
    output:
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta.fai",
    log:
        f"logs/genome-faidx_{config['ref']['build']}_{config['ref']['release']}.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule bwa_index:
    input:
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta",
    output:
        multiext((f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"logs/bwa_index_{config['ref']['build']}_{config['ref']['release']}.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"


rule star_index:
    input:
        fasta=f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta",
        annotation=f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.gtf",
    output:
        directory(f"{assembly_path}star_genome_{config['ref']['build']}_{config['ref']['release']}"),
    threads: 24
    params:
        extra=f"--sjdbGTFfile {assembly_path}{config['ref']['build']}_{config['ref']['release']}.gtf --sjdbOverhang 100",
    log:
        f"logs/star_index_genome_{config['ref']['build']}_{config['ref']['release']}.log",
    cache: True
    wrapper:
        "0.77.0/bio/star/index"

rule rsem_ref:
    input:
        # reference FASTA with either the entire genome or transcript sequences
        f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.gtf",
        reference_genome=f"{assembly_path}{config['ref']['build']}_{config['ref']['release']}.fasta",
    output:
        # one of the index files created and used by RSEM (required)
        seq=f"{assembly_path}rsem_reference_{config['ref']['build']}_{config['ref']['release']}.seq",
        # RSEM produces a number of other files which may optionally be specified as output; these may be provided so that snakemake is aware of them, but the wrapper doesn't do anything with this information other than to verify that the file path prefixes match that of output.seq.
        # for example,
        grp=f"{assembly_path}rsem_reference_{config['ref']['build']}_{config['ref']['release']}.grp",
        ti=f"{assembly_path}rsem_reference_{config['ref']['build']}_{config['ref']['release']}.ti",
    threads: 4
    params:
        # optional additional parameters, for example,
        #extra="--gtf annotations.gtf",
        # if building the index against a reference transcript set
        extra=f"--gtf {assembly_path}{config['ref']['build']}_{config['ref']['release']}.gtf",
        out_ref = f"{assembly_path}rsem_reference_{config['ref']['build']}_{config['ref']['release']}",
    log:
        f"logs/rsem/prepare-reference_{config['ref']['build']}_{config['ref']['release']}.log",
    cache: True
    shell:
        "rsem-prepare-reference --num-threads {threads} {params.extra} {input.reference_genome} {params.out_ref} > {log} 2>&1"