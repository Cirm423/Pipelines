import os

if genecode_assembly:

    rule faToTwoBit_fa:
        input:
            f"{config['resources']}{config['ref']['assembly']}.fa",
        output:
            temp(f"{config['resources']}{config['ref']['assembly']}.2bit"),
        log:
            "logs/browser/fa_to_2bit.log"
        params:
            "" # optional params string
        wrapper:
            "0.78.0/bio/ucsc/faToTwoBit"

    rule twoBitInfo:
        input:
            f"{config['resources']}{config['ref']['assembly']}.2bit"
        output:
            temp(f"{config['resources']}{config['ref']['assembly']}.chrom.sizes.tmp")
        log:
            "logs/browser/chrom.sizes.log"
        params:
            "" # optional params string
        wrapper:
            "0.78.0/bio/ucsc/twoBitInfo"

    rule twoBitInfo_sort:
        input:
            f"{config['resources']}{config['ref']['assembly']}.chrom.sizes.tmp"
        output:
            f"{config['resources']}{config['ref']['assembly']}.chrom.sizes"
        cache: True
        shell:
            "sort -k2rn {input} > {output}"

rule star_str_bw:
    input:
        bam = get_star_bam,
    output:
        temp(path_merged_cond("results/bw_str/?/Signal.Unique.str1.out.bg")),
        temp(path_merged_cond("results/bw_str/?/Signal.Unique.str2.out.bg")),
        temp(path_merged_cond("results/bw_str/?/Signal.UniqueMultiple.str1.out.bg")),
        temp(path_merged_cond("results/bw_str/?/Signal.UniqueMultiple.str2.out.bg")),
    log:
        path_merged_cond("logs/browser/star_str/?.log"),
    params:
    # strandedness
        strand="Stranded",
        path=lambda wc, output: os.path.dirname(output[0]) + "/",
    threads: 8
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode inputAlignmentsFromBAM --inputBAMfile {input.bam} --outWigType bedGraph --outWigStrand {params.strand} --outFileNamePrefix {params.path} --outWigReferencesPrefix chr"

rule BamCoverage_str1:
    input:
        get_star_bam,
        get_star_bam_bai,
    output:
        path_merged_cond("results/browser/?.str1.bw"),
    params:
        norm = config["params"]["bamcoverage"],
        stranded = "" if get_strandedness(units)[0] == 0.5 else "--filterRNAstrand forward"
    log: 
        path_merged_cond("logs/browser/?.BamCoverage.log")
    conda:
        "../envs/deeptools.yaml"
    threads: 12
    shell:
        "bamCoverage -b {input} -o {output} -of bigwig -p {threads} {params.norm} {params.stranded} 2>{log}"

rule BamCoverage_str2:
    input:
        get_star_bam,
        get_star_bam_bai,
    output:
        path_merged_cond("results/browser/?.str2.bw"),
    params:
        norm = config["params"]["bamcoverage"],
        stranded = "--filterRNAstrand reverse"
    log: 
        path_merged_cond("logs/browser/?.BamCoverage.log")
    conda:
        "../envs/deeptools.yaml"
    threads: 12
    shell:
        "bamCoverage -b {input} -o {output} -of bigwig -p {threads} {params.norm} {params.stranded} 2>{log}"



### Old bigwig file code, in case it's useful
# rule star_uns_bw:
#     input:
#         bam = get_star_bam,
#     output:
#         temp(path_merged_cond("results/bw_uns/?/Signal.Unique.str1.out.bg")),
#         temp(path_merged_cond("results/bw_uns/?/Signal.UniqueMultiple.str1.out.bg")),
#     log:
#         path_merged_cond("logs/browser/star_uns/?.log"),
#     params:
#         # strandedness
#         strand="Unstranded",
#         path=lambda wc, output:os.path.dirname(output[0]) + "/",
#     threads: 8
#     conda:
#         "../envs/star.yaml"
#     shell:
#         "STAR --runMode inputAlignmentsFromBAM --inputBAMfile {input.bam} --outWigType bedGraph --outWigStrand {params.strand} --outFileNamePrefix {params.path} --outWigReferencesPrefix chr"

# rule sort_bw_uns:
#     input:
#         path_merged_cond("results/bw_uns/?/Signal.Unique.str1.out.bg"),
#     output:
#         temp(path_merged_cond("results/bw_uns/?/Signal.Unique.str1.out.bg.sorted")),
#     threads: 8
#     conda:
#         "../envs/coreutils.yaml"
#     shell:
#         "sort --parallel={threads} -k1,1 -k2,2n {input} > {output}"

# rule sort_bw_str:
#     input:
#         str1=path_merged_cond("results/bw_str/?/Signal.Unique.str1.out.bg"),
#         str2=path_merged_cond("results/bw_str/?/Signal.Unique.str2.out.bg"),
#     output:
#         str1=temp(path_merged_cond("results/bw_str/?/Signal.Unique.str1.out.bg.sorted")),
#         str2=temp(path_merged_cond("results/bw_str/?/Signal.Unique.str2.out.bg.sorted")),
#     threads: 8
#     conda:
#         "../envs/coreutils.yaml"
#     shell:
#         "sort --parallel={threads} -k1,1 -k2,2n {input.str1} > {output.str1} && sort --parallel={threads} -k1,1 -k2,2n {input.str2} > {output.str2}"

# rule bedGraphToBigWig_str1:
#     input:
#         bedGraph=get_bg_str1,
#         chromsizes=f"{config['resources']}{config['ref']['assembly']}.chrom.sizes",
#     output:
#         path_merged_cond("results/browser/?.str1.bw"),
#     log:
#         path_merged_cond("logs/browser/?.bed-graph_to_big-wig_str1.log"),
#     params:
#         "" # optional params string
#     wrapper:
#         "0.77.0/bio/ucsc/bedGraphToBigWig"

# rule bedGraphToBigWig_str2:
#     input:
#         bedGraph=path_merged_cond("results/bw_str/?/Signal.Unique.str2.out.bg.sorted"),
#         chromsizes=f"{config['resources']}{config['ref']['assembly']}.chrom.sizes",
#     output:
#         path_merged_cond("results/browser/?.str2.bw"),
#     log:
#         path_merged_cond("logs/browser/?.bed-graph_to_big-wig_str1.log"),
#     params:
#         "" # optional params string
#     wrapper:
#         "0.77.0/bio/ucsc/bedGraphToBigWig"