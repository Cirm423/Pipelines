rule fanc_expected:
    input:
        "results/hic/{sample_group}.hic",
    output:
        tsv = "results/matrix_analysis/{sample_group}_expected_contatcs.tsv",
        plot = report("results/matrix_analysis/{sample_group}_distance_decay.pdf",category="Expected values")
    params:
        label = lambda wildcards: f"-l {wildcards.sample_group}",
        extra = config["params"]["fanc"]["analysis"]["expected_params"]
    log:
        "logs/analysis/{sample_group}_expected.log"
    threads: 1
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc expected -p {output.plot} {params.label} -c {params.extra} {input} {output.tsv} 2>{log}"

rule fanc_pca:
    input:
        expand("results/hic/{sample}.hic", sample = samples.index)
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
        hic = "results/hic/{sample_group}.hic",
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

#No insulation based TADs for now
# rule fanc_insulation:
#     input:
#         "results/hic/{sample_group}.hic"
#     output:
#         "results/matrix_analysis/TADs/{sample_group}.insulation",
#     params:
#         extra = config["params"]["fanc"]["analysis"]["TAD_params"]
#     log:
#         "logs/analysis/{sample_group}_TAD.log"
#     threads: 1
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         "fanc insulation {input} {output} {params.extra} 2>{log}"

# rule fanc_insulation_format:
#     input:
#         "results/hic/{sample_group}.hic"
#     output:
#         directory("results/matrix_analysis/TADs/bigwig/{sample_group}"),
#     params:
#         prefix = "results/matrix_analysis/TADs/bigwig/{sample_group}/insulation"
#     log:
#         "logs/analysis/{sample_group}_TAD.log"
#     threads: 1
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         "fanc insulation {input} {params.prefix} -o bigwig 2>{log}"

# rule plot_TAD:
#     input:
#         hic = "results/hic/{sample_group}.hic",
#         insulation = "results/matrix_analysis/TADs/{sample_group}.insulation"
#     output:
#         report("results/matrix_analysis/TADs/{sample_group}.{region}.png",category ="TADs")
#     params:
#         extra = config["params"]["fanc"]["analysis"]["TAD_plot"]
#     log:
#         "logs/analysis/{sample_group}_{region}_compartments_plot.log"
#     threads: 1
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         "fanc plot -o {output} {wildcards.region} -p triangular {input.hic} {params.extra} -p scores {input.insulation} 2>{log}"

rule fanc_loops_annotate:
    input:
        "results/hic/{sample_group}.hic"
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