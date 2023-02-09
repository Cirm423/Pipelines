rule fanc_expected:
    input:
        "results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic",
    output:
        tsv = "results/matrix_analysis/{sample_group}_{chr}.{enzyme}.{fragments}-{resolution}.expected_contatcs.tsv",
        plot = report("results/matrix_analysis/{sample_group}_{chr}.{enzyme}.{fragments}-{resolution}.distance_decay.pdf",category="Expected values")
    params:
        label = lambda wildcards: f"-l {wildcards.sample_group}",
        extra = config["params"]["fanc"]["analysis"]["expected_params"]
    log:
        "logs/analysis/{sample_group}_{chr}.{enzyme}.{fragments}-{resolution}.expected.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc expected -p {output.plot} {params.label} -c {params.extra} {input} {output.tsv} 2>{log}"

rule fanc_pca:
    input:
        expand("results/hic/{sample}.{enzyme}.{fragments}-{resolution}.hic", sample = samples.index, enzyme = enzyme_file, fragments = fragments_file, resolution=analysis_resolution)
    output:
        out = "results/pca/matrix.{enzyme}.{fragments}-{resolution}.pca",
        plot = report("results/pca/matrix.{enzyme}.{fragments}-{resolution}.pca_plot.pdf",category="PCA")
    params:
        label = "-n \"" + '\" \"'.join(samples.index) + "\"",
        extra = config["params"]["fanc"]["analysis"]["pca_params"]
    log:
        "logs/analysis/PCA.{enzyme}.{fragments}-{resolution}.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc pca {params.label} -p {output.plot} {params.extra} {input} {output.out} 2>{log}"

rule fanc_compartments:
    input:
        hic = "results/hic/{group}.{enzyme}.{fragments}-{resolution}.hic",
        fa = f"{assembly_path}{assembly}.fa"
    output:
        AB = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}._compartments.ab",
        eigen = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}.ev.txt",
        domains = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}.domains.bed",
        plot = report("results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}.ab_profile.png", category="AB domains"),
        matrix = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}.enrichment_matrix.txt"
    params:
        extra = config["params"]["fanc"]["analysis"]["AB_params"]
    log:
        "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}_compartments.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc compartments -v {output.eigen} -g {input.fa} -d {output.domains} -e {output.plot} -m {output.matrix} {params.extra} {input.hic} {output.AB} 2>{log}"

rule plot_compartments:
    input:
        AB = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}._compartments.ab",
        eigen = "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}-{resolution}.ev.txt"
    output:
        report("results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}.{region}-{resolution}.compartments.png",category ="AB domains")
    params:
        extra = config["params"]["fanc"]["analysis"]["AB_plot"]
    log:
        "logs/analysis/{group}_{region}.{enzyme}.{fragments}-{resolution}.compartments_plot.log"
    threads: 4
    conda:
        "../envs/fanc.yaml"
    shell:
        "fancplot -o {output} {wildcards.region} -p square {input.AB} {params.extra} -p line {input.eigen} 2>{log}"

# #Only directionality index for TAD calling, not insulation. Fanc directionality does not call TADs!!!
# rule fanc_directionality:
#     input:
#         "results/hic/{group}.{enzyme}.{fragments}-{resolution}.hic"
#     output:
#         "results/matrix_analysis/TADs/{group}.{enzyme}.{fragments}-{resolution}.directionality",
#     params:
#         extra = config["params"]["fanc"]["analysis"]["TAD_params"]
#     log:
#         "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}.TAD.log"
#     threads: 4
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         "fanc insulation {input} {output} {params.extra} 2>{log}"

# rule fanc_directionality_format:
#     input:
#         f"results/matrix_analysis/TADs/{{group}}.{enzyme_file}.{fragments_file}-{analysis_resolution}.directionality"
#     output:
#         directory("results/matrix_analysis/TADs/output/{group}"),
#     params:
#         form = config["params"]["fanc"]["analysis"]["TAD_format"],
#         prefix = lambda wildcards: f"results/matrix_analysis/TADs/output/{wildcards.sample}/{wildcards.sample}.{enzyme_file}.{fragments_file}-{analysis_resolution}.directionality"
#     log:
#         "logs/analysis/{group}_TAD_format.log"
#     threads: 4
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         """
#         mkdir -p {output}
#         fanc insulation {input} {params.prefix} -o {params.form} 2>{log}
#         """

# rule plot_TAD:
#     input:
#         hic = "results/hic/{group}.{enzyme}.{fragments}-{resolution}.hic",
#         directionality = "results/matrix_analysis/TADs/{group}.{enzyme}.{fragments}-{resolution}.directionality"
#     output:
#         report("results/matrix_analysis/TADs/{group}.{enzyme}.{fragments}.{region}-{resolution}.png",category ="TADs")
#     params:
#         extra = config["params"]["fanc"]["analysis"]["TAD_plot"]
#     log:
#         "logs/analysis/{group}_{region}.{enzyme}.{fragments}-{resolution}.directionality_plot.log"
#     threads: 4
#     conda:
#         "../envs/fanc.yaml"
#     shell:
#         "fancplot -o {output} {wildcards.region} -p triangular {input.hic} {params.extra} -p scores {input.directionality} 2>{log}"

rule fanc_loops_annotate:
    input:
        "results/hic/{group}.{enzyme}.{fragments}-{resolution}.hic"
    output:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_annotate"]
    log:
        "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}.loop_annotate.log"
    threads: 24
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} -t {threads} 2>{log}"

rule fanc_loops_filter:
    input:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.loops",
    output:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.filtered_loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_filter"]
    log:
        "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}.loop_filter.log"
    threads: 4
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} -t {threads} 2>{log}"

rule fanc_loops_merge:
    input:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.filtered_loops",
    output:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.merged_loops",
    params:
        extra = config["params"]["fanc"]["analysis"]["Loop_merge"]
    log:
        "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}.loop_merge.log"
    threads: 4
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} {output} {params.extra} -t {threads} 2>{log}"

rule fanc_loops_export:
    input:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.merged_loops",
    output:
        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.merged.bedpe",
    log:
        "logs/analysis/{group}.{enzyme}.{fragments}-{resolution}.loop_export.log"
    threads: 4
    conda:
        "../envs/fanc.yaml"
    shell:
        "fanc loops {input} -b {output} -t {threads} 2>{log}"

rule domaincaller:
    input:
        hic = "results/cooler/{group}.{enzyme}.{fragments}-{resolution}.mcool",
        ini = "results/tadlib/package.done"
    output:
        out = "results/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.tad.bed",
        di_out = "results/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.DIs.bedgraph"
    params:
        extra = config["params"]["fanc"]["analysis"]["domaincaller_extra"],
        uri = lambda wildcards, input: input.hic + "::/resolutions/" + TAD_res
    log:
        "logs/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.log"
    threads: 24
    conda:
        "../envs/tadlib.yaml"
    shell:
        "domaincaller --uri {params.uri} -O {output.out} -D {output.di_out} -p {threads} {params.extra} --logFile {log} --removeCache"

#Add script, check that it works first
rule plot_single_TADs:
    input:
        hic = "results/cooler/{group}.{enzyme}.{fragments}-{resolution}.mcool",
        tad = "results/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.tad.bed",
        loops = "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.merged.bedpe"
        ini = "results/tadlib/package.done"
    output:
        report("results/domaincaller/{group}.{enzyme}.{fragments}.{region}-{resolution}.tad.png",category="TAD calling"),
    params:
        uri = lambda wildcards, input: input.hic + "::/resolutions/" + TAD_res,
        coords = get_coord_params
    log:
        "logs/hitad/{group}.{enzyme}.{fragments}-{resolution}.{region}_TAD-plot.log"     
    threads: 1
    conda:
        "../envs/tadlib.yaml"
    script:
        "../scripts/plot_TADs.py"

rule create_hitad_meta:
    input:
        files = get_tadlib_input,
        ini = "results/tadlib/package.done"
    output:
        "results/hitad/{group}.{enzyme}.{fragments}-{resolution}.hitad_meta.txt"
    params:
        TAD_res
    log:
        "logs/hitad/{group}.{enzyme}.{fragments}-{resolution}_meta.log"
    run:
        with open(output[0],"w") as fout:
            fout.write(f"res:{params[0]}\n")
            for fin in input.files:
                label = fin.split("/")[-1].split(f".{wildcards.enzyme}")[0]
                fout.write(f"  {label}:{fin}::/resolutions/{TAD_res}\n")

rule hitad:
    input:
        meta = "results/hitad/{group}.{enzyme}.{fragments}-{resolution}.hitad_meta.txt",
        ini = "results/tadlib/package.done"
    output:
        "results/hitad/{group}.{enzyme}.{fragments}-{resolution}.hitad.txt"
    params:
        config["params"]["fanc"]["analysis"]["hitad_extra"]
    log:
        "logs/hitad/{group}.{enzyme}.{fragments}-{resolution}_hitad.log"
    threads: 24
    conda:
        "../envs/tadlib.yaml"
    shell:
        "hitad -O {output} -d {input.meta} --logFile {log} {params} --removeCache -p {threads}"

rule plot_hierarchical_TADs:
    input:
        hitad = get_hitad_output,
        hic = "results/cooler/{sample}.{enzyme}.{fragments}-{resolution}.mcool",
        ini = "results/tadlib/package.done"
    output:
        report("results/hitad/{sample}.{enzyme}.{fragments}.{region}-{resolution}.hitad.png",category="TAD calling")
    params:
        uri = lambda wildcards, input: input.hic + "::/resolutions/" + TAD_res,
        extra = config["params"]["fanc"]["analysis"]["TAD_plot"],
        coords = get_coord_params
    log:
        "logs/hitad/{sample}.{enzyme}.{fragments}-{resolution}.{region}_hitad-plot.log"
    conda:
        "../envs/tadlib.yaml"
    shell:
        "tad-plot -p {params.uri} -T {input.hitad} -O {output} {params.coords} {params.extra} 2>{log}"

rule multi_TAD_browser:
    input:
        hitad = get_hitad_output,
        hic = "results/cooler/{sample}.{enzyme}.{fragments}-{resolution}.mcool",
        ini = "results/tadlib/package.done"
    output:
        "results/hitad/{sample}.{enzyme}.{fragments}-{resolution}.hitad.DIs.BedGraph"
    params:
        uri = lambda wildcards, input: input.hic + "::/resolutions/" + TAD_res
    log:
        "logs/hitad/{sample}.{enzyme}.{fragments}-{resolution}_DI.log"
    conda:
        "../envs/tadlib.yaml"
    shell:
        "output-DI -p {params.uri} -O {output} 2>{log}"