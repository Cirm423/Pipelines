import os

rule bamtobed_TE_single:
    input:
        get_star_bam,
    output:
        "results/TE_single/{sample}/Aligned.out.bed",
    log:
        "logs/TE_single/bamtobed/{sample}.log",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -split -i {input} > {output} 2> {log}"

rule TE_single_signal:
    input:
        bed="results/TE_single/{sample}/Aligned.out.bed",
        fac="results/Factor",
        annot=f"{config['resources']}{config['ref']['assembly']}.rmsk.bed",
        log=get_star_log,
    output:
        tmp=temp("results/TE_single/{sample}/Aligned.out.bed.sorted"),
        out="results/TE_single/{sample}/uniq-TE-individual.signal",
    log:
        "logs/TE_single/signal/{sample}.log",
    threads: 4
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bash workflow/scripts/single/TE_single_signal.bash {input.bed} {input.fac} {input.annot} {input.log} {output.tmp} {output.out} 2 > {log}"

rule filter_single_TE:
    input:
        "results/TE_single/{sample}/uniq-TE-individual.signal",
    output:
        "results/TE_single/{sample}/uniq-TE-individual.signal.filtered",
    params:
        filt=config["TE_single"]["filter"],
    shell:
        r"""awk -F"\t" 'BEGIN{{OFS="\t"}} {{if ($7>{params.filt}) print }}' {input} > {output}"""

rule TE_single_diff:
    input:
        treat = expand("results/TE_single/{sample}/uniq-TE-individual.signal.filtered",sample=samples.loc[samples["condition"]!="control"]["sample_name"]),
        control = expand("results/TE_single/{sample}/uniq-TE-individual.signal.filtered",sample=samples.loc[samples["condition"]=="control"]["sample_name"]),
    output:
        expand("results/TE_single/{sample}_vs_{control}.tsv",sample=samples.loc[samples["condition"]!="control"]["sample_name"], control = samples.loc[samples["condition"]=="control"]["sample_name"]),
        temp("results/TE_single/counts.temp"),
        temp("results/TE_single/counts_diff.temp"),
    run:
        for file in input.treat:
            file_name = os.path.basename(os.path.dirname(file))
            control_name = os.path.basename(os.path.dirname(input.control[0]))
            shell(r"""paste {file} {input.control}  | awk 'BEGIN{{OFS="\t"}} {{if($15>0) {{print $1,$2,$3,$4,$5,$6,$7,$15,$8,$16,$7/$15}} else {{print $1,$2,$3,$4,$5,$6,$7,"1",$8,$16,$7/1}}}}' > results/TE_single/counts.temp""")
            shell("""
            source activate R
            Rscript workflow/scripts/pois-test_uniq-TE-individual.R
            source deactivate
            """)
            shell(r"""echo -e "Chr\tStart\tEnd\tRegion_ID\tMappability\tStrand\t{file_name}_NormalizedCounts\t{control_name}_NormalizedCounts\t{file_name}_rpkm\t{control_name}_rpkm\tFoldChange\tp-value\tFDR" > results/TE_single/{file_name}_vs_{control_name}.tsv""")
            shell("cat results/TE_single/counts_diff.temp >> results/TE_single/{file_name}_vs_{control_name}.tsv")