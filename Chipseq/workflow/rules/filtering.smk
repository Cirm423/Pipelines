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