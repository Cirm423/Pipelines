rule methylkit:
    input:
        meth = get_methylkit_input,
        annot = f"{assembly_path}{assembly}.annotation.bed"
    output:
        temp(directory("results/MethylDB")),
        RData = "results/diff_meth/methykit.RData",
        CpG_methylation = report("results/diff_meth/plots/CpG_methylation_percent.pdf", category="Differential Methylation"),
        CpG_coverage = report("results/diff_meth/plots/CpG_coverage.pdf", category="Differential Methylation"),
        correlation = report("results/diff_meth/plots/Sample_correlation.pdf", category="Differential Methylation"),
        cluster = report("results/diff_meth/plots/Sample_clustering.pdf", category="Differential Methylation"),
        screen = report("results/diff_meth/plots/PCA_screen.pdf", category="Differential Methylation"),
        PCA = report("results/diff_meth/plots/PCA.pdf", category="Differential Methylation"),
        hyper = "results/diff_meth/CpG_hypermethylated_25p.tsv",
        hypo = "results/diff_meth/CpG_hypomethylated_25p.tsv",
        all_diff = "results/diff_meth/CpG_all_methylated_diff_25p.tsv",
        chr_diff = "results/diff_meth/CpG_methylated_by_chr_25p.tsv",
        annotation = "results/diff_meth/CpG_methylated_annotation_25p.tsv"
    params:
        mode = config["params"]["mode"],
        assembly = assembly,
        treatment = [1 if x == "treatment" else 0 for x in samples["group"]],
        samples = config["samples"],
        #Changes the file input to methylkit
        mincov = config["params"]["methylkit"]["mincov"],
        min_group = config["params"]["methylkit"]["min_group"]
    threads: 24
    conda:
        "../envs/methylkit.yaml"
    script:
        "../scripts/methylkit.R"