rule deseq2_init:
    input:
        counts=expand("results/rsem/{star_lib}/{sample}/mapped.genes.results",star_lib=star_lib,sample=samples.sample_name),
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["params"]["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/plots/pca.svg", caption = "../report/pca.rst", category = "PCA"),
    params:
        pca_labels=config["params"]["pca"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log",
    script:
        "../scripts/plot-pca.R"


rule deseq2:
    input:
        "results/deseq2/all.rds",
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        ma_plot=report("results/plots/{contrast}.ma-plot.svg", caption = "../report/ma.rst", category = "MA plots"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule deseq2_convert:
    input:
        "results/diffexp/{contrast}.diffexp.tsv",
    output:
        report(
            "results/diffexp/{contrast}.diffexp.gene_symbol.tsv", caption = "../report/diffexp.rst", category = "Differential Expression tables"
        ),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2_ID_conv.yaml"
    log:
        "logs/deseq2/{contrast}.ID_conversion.log",
    threads: 2
    script:
        "../scripts/DESeq2_convertID.R"

rule enhanced_volcano:
    input:
        "results/diffexp/{contrast}.diffexp.gene_symbol.tsv",
    output:
        report("results/plots/{contrast}.volcano_plot.svg", category = "Volcano plots"),
    threads: 1
    log:
        "logs/enhanced-volcano/{contrast}.tsv.log",
    params:
        contrast=get_contrast,
        #extra="lab='gene', x='log2FoldChange', y='padj'",
        #width=1024,  # Optional PNG width
        #height=768,  # Optional PNG height
    conda:
        "../envs/enhanced_volcano.yaml"
    script:
        "../scripts/enhanced_volcano.R"

rule deseq2_init_TE:
    input:
        bam=expand("results/filtered/{star_lib}/{sample}.filtered.sortedByCoord.out.bam",star_lib=star_lib,sample=samples.sample_name),
        gtf=f"{assembly_path}{assembly}.rmsk.gtf",
    output:
        "results/deseq2/TE_all.rds",
        "results/deseq2/TE_normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["params"]["diffexp"]["model"],
        end=get_deseq2_end,
        filt=config["params"]["diffexp"]["TE"]["filter"],
    conda:
        "../envs/deseq2_TE.yaml"
    log:
        "logs/deseq2/init_TE.log",
    threads: 24
    script:
        "../scripts/DESeq2-TE.R"

rule pca_TE:
    input:
        "results/deseq2/TE_all.rds",
    output:
        report("results/plots/TE_pca.svg", caption = "../report/pca_TE.rst", category = "PCA"),
    params:
        pca_labels=config["params"]["pca"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca_TE.log",
    script:
        "../scripts/plot-pca.R"

rule deseq2_TE:
    input:
        "results/deseq2/TE_all.rds",
    output:
        table="results/diffexp/{contrast}.diffexp.TE.tsv",
        ma_plot=report("results/plots/{contrast}.ma-plot.TE.svg", caption = "../report/ma_TE.rst", category = "MA plots"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2_TE.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.TE.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule deseq2_convert_TE:
    input:
        "results/diffexp/{contrast}.diffexp.TE.tsv",
    output:
        report(
            "results/diffexp/{contrast}.diffexp.gene_symbol.TE.tsv", caption = "../report/diffexp_TE.rst", category = "Differential expression tables"
        ),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2_ID_conv.yaml"
    log:
        "logs/deseq2/{contrast}.ID_conversion.TE.log",
    threads: 2
    script:
        "../scripts/DESeq2_convertID.R"

rule enhanced_volcano_TE:
    input:
        "results/diffexp/{contrast}.diffexp.gene_symbol.TE.tsv",
    output:
        report("results/plots/{contrast}.volcano_plot.TE.svg", category = "Volcano plots"),
    threads: 1
    log:
        "logs/enhanced-volcano/{contrast}.tsv.TE.log",
    params:
        contrast=get_contrast,
        #extra="lab='gene', x='log2FoldChange', y='padj'",
        #width=1024,  # Optional PNG width
        #height=768,  # Optional PNG height
    conda:
        "../envs/enhanced_volcano.yaml"
    script:
        "../scripts/enhanced_volcano.R"
