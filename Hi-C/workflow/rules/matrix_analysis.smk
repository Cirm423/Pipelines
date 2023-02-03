rule fanc_expected:
    input:
        f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic",
    output:
        tsv = "results/matrix_analysis/{sample_group}_{chr}_expected_contatcs.tsv",
        plot = report("results/matrix_analysis/{sample_group}_{chr}_distance_decay.pdf",category="Expected values")
    params:
        label = lambda wildcards: f"-l {wildcards.sample_group}",
        extra = config["params"]["fanc"]["analysis"]["expected_params"]
    log:
        "logs/analysis/{sample_group}_{chr}_expected.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc expected -p {output.plot} {params.label} -c {params.extra} {input} {output.tsv} 2>{log}"

rule fanc_pca:
    input:
        expand(f"results/hic/{{sample}}.{enzyme_file}.{fragments_file}.hic", sample = samples.index)
    output:
        out = "results/pca/matrix.pca",
        plot = report("results/pca/matrix_pca_plot.pdf",category="PCA")
    params:
        label = "-n \"" + '\" \"'.join(samples.index) + "\"",
        extra = config["params"]["fanc"]["analysis"]["pca_params"]
    log:
        "logs/analysis/PCA.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc pca {params.label} -p {output.plot} {params.extra} {input} {output.out} 2>{log}"

rule fanc_compartments:
    input:
        hic = f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic",
        fa = f"{assembly_path}{assembly}.fa"
    output:
        AB = "results/matrix_analysis/compartments/{sample_group}_compartments.ab",
        eigen = "results/matrix_analysis/compartments/{sample_group}.ev.txt",
        domains = "results/matrix_analysis/compartments/{sample_group}_domains.bed",
        plot = report("results/matrix_analysis/compartments/{sample_group}_ab_profile.png", category="AB domains"),
        matrix = "results/matrix_analysis/compartments/{sample_group}_enrichment_matrix.txt"
    params:
        extra = config["params"]["fanc"]["analysis"]["AB_params"]
    log:
        "logs/analysis/{sample_group}_compartments.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc compartments -v {output.eigen} -g {input.fa} -d {output.domains} -e {output.plot} -m {output.matrix} {params.extra} {input.hic} {output.AB} 2>{log}"

rule plot_compartments:
    input:
        AB = "results/matrix_analysis/compartments/{sample_group}_compartments.ab",
        eigen = "results/matrix_analysis/compartments/{sample_group}.ev.txt"
    output:
        report("results/matrix_analysis/compartments/{sample_group}.{region}.png",category ="AB domains")
    params:
        extra = config["params"]["fanc"]["analysis"]["AB_plot"]
    log:
        "logs/analysis/{sample_group}_{region}_compartments_plot.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fancplot -o {output} {wildcards.region} -p square {input.AB} {params.extra} -p line {input.eigen} 2>{log}"

#Only directionality index for TAD calling, not insulation.
rule fanc_directionality:
    input:
        f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic"
    output:
        "results/matrix_analysis/TADs/{sample_group}.directionality",
    params:
        extra = config["params"]["fanc"]["analysis"]["TAD_params"]
    log:
        "logs/analysis/{sample_group}_TAD.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc insulation {input} {output} {params.extra} 2>{log}"

rule fanc_directionality_format:
    input:
        "results/matrix_analysis/TADs/{sample_group}.directionality"
    output:
        directory("results/matrix_analysis/TADs/output/{sample_group}"),
    params:
        form = config["params"]["fanc"]["analysis"]["TAD_format"],
        prefix = lambda wildcards: f"results/matrix_analysis/TADs/output/{wildcards.sample_group}/{wildcards.sample_group}.directionality"
    log:
        "logs/analysis/{sample_group}_TAD_format.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc insulation {input} {params.prefix} -o {params.form} 2>{log}"

rule plot_TAD:
    input:
        hic = f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic",
        directionality = "results/matrix_analysis/TADs/{sample_group}.directionality"
    output:
        report("results/matrix_analysis/TADs/{sample_group}.{region}.png",category ="TADs")
    params:
        extra = config["params"]["fanc"]["analysis"]["TAD_plot"]
    log:
        "logs/analysis/{sample_group}_{region}_directionality_plot.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc plot -o {output} {wildcards.region} -p triangular {input.hic} {params.extra} -p scores {input.directionality} 2>{log}"

rule fanc_loops_annotate:
    input:
        f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic"
    output:
        "results/matrix_analysis/loops/{sample_group}.loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_annotate"]
    log:
        "logs/analysis/{sample_group}_loop_annotate.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} 2>{log}"

rule fanc_loops_filter:
    input:
        "results/matrix_analysis/loops/{sample_group}.loops",
    output:
        "results/matrix_analysis/loops/{sample_group}_filtered.loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_filter"]
    log:
        "logs/analysis/{sample_group}_loop_filter.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} 2>{log}"

rule fanc_loops_merge:
    input:
        "results/matrix_analysis/loops/{sample_group}_filtered.loops",
    output:
        "results/matrix_analysis/loops/{sample_group}_merged.loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_merge"]
    log:
        "logs/analysis/{sample_group}_loop_merge.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} 2>{log}"

rule fanc_loops_export:
    input:
        "results/matrix_analysis/loops/{sample_group}_merged.loops",
    output:
        "results/matrix_analysis/loops/{sample_group}_merged.bedpe",
    log:
        "logs/analysis/{sample_group}_loop_export.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} -b {output} 2>{log}"

#Domaincaller replaced by the fanc function above

# rule domaincaller:
#     input:
#         hic = "results/cooler/{sample_group}.cooler.mcool",
#         ini = "results/domaincaller/package.done"
#     output:
#         out = "results/domaincaller/{sample_group}.cooler.domains",
#         di_out = "results/domaincaller/{sample_group}.cooler.di_domains"
#     params:
#         extra = config["params"]["fanc"]["analysis"]["domaincaller_extra"],
#         uri = lambda wildcards, input: input.hic + "::/resolutions/" + config["params"]["fanc"]["analysis"]["domaincaller_res"]
#     log:
#         "logs/domaincaller/{sample_group}.log"
#     threads: 24
#     conda:
#         "../envs/domaincaller.yaml"
#     shell:
#         "domaincaller --uri {params.uri} -O {output.out} -D {output.di_out} -p {threads} {params.extra} --logFile {log}"