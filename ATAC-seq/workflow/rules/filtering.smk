rule samtools_view_filter:
    input:
        "results/merged/{sample}.bam", 
        f"{config['resources']['path']}{config['resources']['ref']['assembly']}.blacklist.sorted.complement"
    output:
        temp("results/sam-view/{sample}.bam")
    params:
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # the -L option is automatically activated if a path to a blacklist of the given genome exists in the
        # downloaded "resources/ref/igenomes.yaml" or has been provided via the "config/config.yaml"
        # parameter "config['resources']['ref']['blacklist']"
        extra=lambda wc, input: "{params} {blacklist}".format(
            params=config["params"]["samtools-view-se"] if config["single_end"] else config["params"]["samtools-view-pe"],
            blacklist="" if len(input) == 1 else "-L {}".format(list(input)[1])
        )
    log:
        "logs/samtools-view/{sample}.log"
    wrapper:
        "v1.3.1/bio/samtools/view"

rule bamtools_filter_json:
    input:
        "results/sam-view/{sample}.bam"
    output:
        temp("results/bamtools_filtered/{sample}.bam")
    params:
          # filters mismatches in all reads and filters pe-reads within a size range given in json-file
        json="config/{}_bamtools_filtering_rules.json".format("se" if config["single_end"] else "pe"),
        region=""
    log:
        "logs/filtered/{sample}.log"
    wrapper:
        "v1.3.1/bio/bamtools/filter_json"

rule samtools_sort:
    input:
        "results/bamtools_filtered/{sample}.bam"
    output:
        temp("results/bamtools_filtered/{sample}.sorted.bam")
    params:
        extra=""
    log:
        "logs/bamtools_filtered/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

#TODO for later: customize and substitute rm_orphan_pe_bam.py with some existing tool
rule orphan_remove:
    input:
        "results/bamtools_filtered/{sample}.sorted.bam"
    output:
        bam=temp("results/orph_rm_pe/{sample}.bam"),
        qc="results/orph_rm_pe/{sample}_bampe_rm_orphan.log"
    params:
        "--only_fr_pairs"
    log:
        "logs/orph_rm_pe/{sample}.log"
    conda:
        "../envs/pysam.yaml"
    shell:
        " workflow/scripts/rm_orphan_pe_bam.py {input} {output.bam} {params} 2> {log}"

rule samtools_sort_pe:
    input:
         "results/orph_rm_pe/{sample}.bam"
    output:
        temp("results/orph_rm_pe/{sample}.sorted.bam")
    params:
        extra=""
    log:
        "logs/orph_rm_pe/{sample}.sorted.log"
    threads:  # Samtools takes additional threads through its option -@
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

#Instead of linking the files like in ChIP, here they're sorted to queryname for genrich
rule merge_se_pe:
    input:
        get_se_pe_branches_input
    output:
        "results/filtered/{sample}.sorted.bam"
    params:
        ""
    log:
        "logs/filtered/{sample}.sorted.log"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -sr {input} {output}"