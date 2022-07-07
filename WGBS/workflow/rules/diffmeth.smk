rule methylkit:
    input:
        meth = get_methylkit_input,
        annot = f"{assembly_path}{assembly}.annotation.bed"
    output:
        temp(directory("results/MethylDB")),
        CpG_methylation = report("results/diff_meth/CpG_methylation_percent.pdf", category="Differential Methylation"),
        CpG_coverage = report("results/diff_meth/CpG_coverage.pdf", category="Differential Methylation"),
        correlation = report("results/diff_meth/Sample_correlation.pdf", category="Differential Methylation"),
        cluster = report("results/diff_meth/Sample_clustering.pdf", category="Differential Methylation"),
        screen = report("results/diff_meth/PCA_screen.pdf", category="Differential Methylation"),
        PCA = report("results/diff_meth/PCA.pdf", category="Differential Methylation"),
    params:
        mode = config["params"]["mode"],
        assembly = assembly,
        treatment = [1 if x == "treatment" else 0 for x in samples["group"]],
        samples = config["samples"],
        #Changes the file input to methylkit
        mincov = config["params"]["methylkit"]["mincov"],
        min_group = config["params"]["methylkit"]["min_group"]
    threads: 12
    conda:
        "../envs/methylkit.yaml"
    script:
        "../scripts/methylkit.R"