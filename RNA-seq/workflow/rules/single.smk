rule get_factor:
    input:
        expand("results/star/{star_lib}/{sample}/Log.final.out",star_lib=star_lib,sample=samples.sample_name)
    output:
        "results/Factor"
    script:
        "../scripts/single/find-factor.py"

rule gene_name:
    input:
        f"{assembly_path}{assembly}.annotation.gtf",
    output:
        f"{assembly_path}{assembly}.gtf.gene_name",
    cache: True
    log:
        "logs/gene_name.log",
    shell:
        "bash workflow/scripts/single/gene_name.bash {input} {output} 2 > {log}"

rule prepare_single:
    input:
        treat = expand("results/rsem/{star_lib}/{{sample}}/mapped.genes.results",star_lib=star_lib),
        id_name = f"{assembly_path}{assembly}.gtf.gene_name",
        factor = "results/Factor",
    output:
        "results/single/{sample}.signal"
    shell:
        "bash workflow/scripts/single/prepare_single.bash {input.treat} {input.id_name} {input.factor} {output}"

rule counts_single:
    input:
        treat = expand("results/single/{sample}.signal",sample=samples.loc[samples["condition"]!="control"]["sample_name"]),
        control = expand("results/single/{sample}.signal",sample=samples.loc[samples["condition"]=="control"]["sample_name"]),
    output:
        expand("results/single/{sample}_vs_{control}.tsv",sample=samples.loc[samples["condition"]!="control"]["sample_name"], control = samples.loc[samples["condition"]=="control"]["sample_name"]),
        temp("results/single/counts.temp"),
        temp("results/single/counts_diff.temp"),
    params:
        config["params"]["single"]["alternative"]
    log:
        "logs/single/counts_single.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/single/counts_single.py"