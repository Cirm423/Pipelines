rule deseq2_init:
    input:
        counts=lambda wc: get_star_output_all_units(wc, fi='rsem'),
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
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
        report("results/pca.svg", caption = "../report/pca.rst", category = "PCA"),
    params:
        pca_labels=config["pca"]["labels"],
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
        table=report(
            "results/diffexp/{contrast}.diffexp.tsv", caption = "../report/diffexp.rst", category = "Differential expression"
        ),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", caption = "../report/ma.rst", category = "Differential expression"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"


rule deseq2_init_TE:
    input:
        bam=lambda wc: get_star_output_all_units(wc, fi='bam'),
        gtf=f"{assembly_path}{assembly}.rmsk.gtf",
    output:
        "results/deseq2/TE_all.rds",
        "results/deseq2/TE_normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
        end=get_deseq2_end,
        filt=config["diffexp"]["TE"]["filter"],
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
        report("results/TE_pca.svg", caption = "../report/pca_TE.rst", category = "PCA"),
    params:
        pca_labels=config["pca"]["labels"],
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
        table=report(
            "results/diffexp/{contrast}.diffexp.TE.tsv", caption = "../report/diffexp_TE.rst", category = "Differential expression"
        ),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.TE.svg", caption = "../report/ma_TE.rst", category = "Differential expression"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2_TE.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.TE.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"
