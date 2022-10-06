rule samtools_view_filter:
    input:
        get_samtools_view_filter_input
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
        "results/bamtools_filtered/{sample}.sorted.bam"
    params:
        extra=""
    log:
        "logs/bamtools_filtered/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "v1.3.1/bio/samtools/sort"

rule bamtobed:
    input:
        "results/bamtools_filtered/{sample}.sorted.bam",
    output:
        temp("results/bamtools_filtered/{sample}.bed"),
    log:
        "logs/bamtobed/{sample}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -bedpe -i {input} > {output} 2> {log}"

rule clean_bed:
    input:
        "results/bamtools_filtered/{sample}.bed",
    output:
        temp("results/bamtools_filtered/{sample}_clean.bed"),
    log:
        "logs/bamtobed/{sample}.clean.log"
    shell:
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input} > {output} 2>{log}"

rule fragment_bed:
    input:
        "results/bamtools_filtered/{sample}_clean.bed",
    output:
        temp("results/bamtools_filtered/{sample}_fragments.bed"),
    log:
        "logs/bamtobed/{sample}.fragment.log"
    shell:
        "cut -f 1,2,6 {input} | sort -k1,1 -k2,2n -k3,3n > {output}"

rule genomecov_bed_spike:
    input:
        bed="results/bamtools_filtered/{sample}_fragments.bed",
        ref=f"{assembly_path}{assembly}.chrom.sizes",
        flag_stats="results/mapped/{sample}_spike-in.bam.flagstat",
    output:
        "results/bed_graph/{sample}_normalized.bedgraph"
    log:
        "logs/bed_graph/{sample}.log"
    params:
        lambda w, input:
            "-bg -scale $(grep -m 1 'mapped (' {flagstats_file} | awk '{{print 10000/$1}}')".format(
            flagstats_file=input.flag_stats,
        )
    wrapper:
        "v1.3.1/bio/bedtools/genomecov"

rule genomecov_bed_no_spike:
    input:
        bed="results/bamtools_filtered/{sample}_fragments.bed",
        ref=f"{assembly_path}{assembly}.chrom.sizes"
    output:
        "results/bed_graph/{sample}.bedgraph"
    log:
        "logs/bed_graph/{sample}.log"
    params:
        "-bg"
    wrapper:
        "v1.3.1/bio/bedtools/genomecov"