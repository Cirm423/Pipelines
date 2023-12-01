import os

rule bamtobed_TE_single:
    input:
        expand("results/star/{lib}/{sample}/Aligned.sortedByCoord.out.bam",lib=star_lib,sample=samples.sample_name),
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
        annot=f"{assembly_path}{assembly}.rmsk.bed",
        log=expand("results/star/{star_lib}/{{sample}}/Log.final.out",star_lib=star_lib),
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
        filt=config["params"]["TE_single"]["filter"],
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
    log:
        "logs/TE_single/TE_diff.log"
    params:
        config["params"]["TE_single"]["alternative"]
    run:
        for file in input.treat:
            file_name = os.path.basename(os.path.dirname(file))
            control_name = os.path.basename(os.path.dirname(input.control[0]))
            shell(r"""paste {file} {input.control}  | awk 'BEGIN{{OFS="\t"}} {{if($15>0) {{print $1,$2,$3,$4,$5,$6,$7,$15,$8,$16,$7/$15}} else {{print $1,$2,$3,$4,$5,$6,$7,"1",$8,$16,$7/1}}}}' > results/TE_single/counts.temp 2>{log}""")
            shell("Rscript workflow/scripts/pois-test_uniq-TE-individual.R {params}")
            shell(r"""echo -e "Chr\tStart\tEnd\tRegion_ID\tMappability\tStrand\t{file_name}_NormalizedCounts\t{control_name}_NormalizedCounts\t{file_name}_rpkm\t{control_name}_rpkm\tFoldChange\tp-value\tFDR" > results/TE_single/{file_name}_vs_{control_name}.tsv 2>>{log}""")
            shell("cat results/TE_single/counts_diff.temp >> results/TE_single/{file_name}_vs_{control_name}.tsv 2>>{log}")