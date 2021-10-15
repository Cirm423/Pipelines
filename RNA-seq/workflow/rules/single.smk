import os

rule get_factor:
    input:
        lambda wc: get_star_output_all_units(wc, fi="log"),
    output:
        "results/Factor"
    script:
        "../scripts/single/find-factor.py"

rule gene_name:
    input:
        f"{config['resources']}{config['ref']['assembly']}.annotation.gtf",
    output:
        f"{config['resources']}{config['ref']['assembly']}.gtf.gene_name",
    cache: True
    log:
        "logs/gene_name.log",
    shell:
        "workflow/scripts/single/gene_name.bash {input} {output} 2 > {log}"

rule prepare_single:
    input:
        treat = get_single_input,
        id_name = f"{config['resources']}{config['ref']['assembly']}.gtf.gene_name",
        factor = "results/Factor",
    output:
        "results/single/{sample}.signal"
    shell:
        "workflow/scripts/single/prepare_single.bash {input.treat} {input.id_name} {input.factor} {output}"

rule counts_single:
    input:
        treat = expand("results/single/{sample}.signal",sample=samples.loc[samples["condition"]!="control"]["sample_name"]),
        control = expand("results/single/{sample}.signal",sample=samples.loc[samples["condition"]=="control"]["sample_name"]),
    output:
        expand("results/single/{sample}_vs_{control}.tsv",sample=samples.loc[samples["condition"]!="control"]["sample_name"], control = samples.loc[samples["condition"]=="control"]["sample_name"]),
        temp("results/single/counts.temp"),
        temp("results/single/counts_diff.temp"),
    run:
        for file in input.treat:
            file_name = os.path.basename(file).split(".")[0]
            control_name = os.path.basename(input.control[0]).split(".")[0]
            shell("""paste {file} {input.control}  | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$12,$6,$13,$7,$14,$5/$12}}' > results/single/counts.temp""")
            shell("""
            source activate R
            Rscript workflow/scripts/pois-test_RSEM.R
            source deactivate
            """)
            shell("""echo -e "GeneName\tGeneID\tlength\teffective_length\t{file_name}_NormalizedCounts\t{control_name}_NormalizedCounts\t{file_name}_TPM\t{control_name}_TPM\t{file_name}_FPKM\t{control_name}_FPKM\tFoldChange\tp-value\tFDR" > results/single/{file_name}_vs_{control_name}.tsv""")
            shell("cat results/single/counts_diff.temp >> results/single/{file_name}_vs_{control_name}.tsv")