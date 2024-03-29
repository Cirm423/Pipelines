from snakemake.utils import validate
import pandas as pd
import os
from smart_open import open
import yaml

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype = str).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

build = config["resources"]["ref"]["assembly"]
chromosome = config["resources"]["ref"]["chromosome"]

assembly = config["resources"]["ref"]["assembly"]
assembly_path = config['resources']['path'] + config['resources']['ref']['assembly'] + "/"

#Check that the control samples are actually in the sample table

assert all(samples[pd.notnull(samples["control"])]["control"].isin(samples.index)), "One or more of the control samples are missing"

##### wildcard constraints #####

wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])

####### helpers ###########

def is_single_end(sample, unit):
    """Determine whether unit is single-end."""
    fq2_present = pd.isnull(units.loc[(sample, unit), "fq2"])
    if isinstance(fq2_present, pd.core.series.Series):
        # if this is the case, get_fastqs cannot work properly
        raise ValueError(
            f"Multiple fq2 entries found for sample-unit combination {sample}-{unit}.\n"
            "This is most likely due to a faulty units.tsv file, e.g. "
            "a unit name is used twice for the same sample.\n"
            "Try checking your units.tsv for duplicates."
        )
    return fq2_present

def has_only_sra_accession(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq1"]) and pd.isnull(units.loc[(sample, unit), "fq2"]) \
           and not pd.isnull(units.loc[(sample, unit), "sra_accession"])

def is_sra_se(sample, unit):
    return has_only_sra_accession(sample, unit) and config["single_end"]

def is_sra_pe(sample, unit):
    return has_only_sra_accession(sample, unit) and not config["single_end"]

def get_individual_fastq(wildcards):
    """Get individual raw FASTQ files from unit sheet, based on a read (end) wildcard"""
    if ( wildcards.read == "0" or wildcards.read == "1" ):
        if is_sra_se(wildcards.sample, wildcards.unit):
            return expand("resources/sra-se-reads/{accession}.fastq.gz",
                              accession=units.loc[ (wildcards.sample, wildcards.unit), "sra_accession" ])
        elif is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("resources/sra-pe-reads/{accession}_1.fastq.gz",
                              accession=units.loc[ (wildcards.sample, wildcards.unit), "sra_accession" ])
        else:
            return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    elif wildcards.read == "2":
        if is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("resources/sra-pe-reads/{accession}_2.fastq.gz",
                          accession=units.loc[ (wildcards.sample, wildcards.unit), "sra_accession" ])
        else:
            return units.loc[ (wildcards.sample, wildcards.unit), "fq2" ]

def get_individual_trimmed_fastq(wildcards):
    if wildcards.read == "0":
        return expand("results/trimmed/{sample}-{unit}.fastq.gz",
                    sample = wildcards.sample, unit = wildcards.unit)
    elif wildcards.read == "1":
        return expand("results/trimmed/{sample}-{unit}_1.fastq.gz",
                        sample = wildcards.sample, unit = wildcards.unit)
    elif wildcards.read == "2":
        return expand("results/trimmed/{sample}-{unit}_2.fastq.gz",
                        sample = wildcards.sample, unit = wildcards.unit)

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_sra_se(wildcards.sample, wildcards.unit):
        return expand("resources/sra-se-reads/{accession}.fastq.gz",
                          accession=units.loc[ (wildcards.sample, wildcards.unit), "sra_accession" ])
    elif is_sra_pe(wildcards.sample, wildcards.unit):
        return expand(["resources/sra-pe-reads/{accession}_1.fastq.gz", "resources/sra-pe-reads/{accession}_2.fastq.gz"],
                          accession=units.loc[ (wildcards.sample, wildcards.unit), "sra_accession" ])
    elif is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    else:
        u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
        return [ f"{u.fq1}", f"{u.fq2}" ]

def is_control(sample):
    control = samples.loc[sample]["control"]
    return pd.isna(control) or pd.isnull(control)

def get_sample_control_peak_combinations_list():
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample):
            sam_contr.extend(expand(["{sample}-{control}.{peak}"], sample = sample, control = samples.loc[sample]["control"], peak = config["params"]["peak-analysis"]))
    return sam_contr

def get_peaks_count_plot_input():
    return expand(
        "results/macs2_callpeak/peaks_count/{sam_contr_peak}.peaks_count.tsv",
        sam_contr_peak = get_sample_control_peak_combinations_list()
    )

def get_frip_score_input():
    return expand(
        "results/bedtools_intersect/{sam_contr_peak}.peaks_frip.tsv",
        sam_contr_peak = get_sample_control_peak_combinations_list()
    )

def get_macs2_peaks():
    return expand(
        "results/macs2_callpeak/{sam_contr_peak}_peaks.{peak}Peak",
        sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
    )

def get_sample_control_peak_combinations_list_ab(antibody):
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample) and samples.loc[sample]["antibody"]==antibody:
            sam_contr.extend(expand(["{sample}-{control}.{peak}"], sample = sample, control = samples.loc[sample]["control"], peak = config["params"]["peak-analysis"]))
    return sam_contr

def get_sample_control_peak_combinations_list_ab_nop(antibody):
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample) and samples.loc[sample]["antibody"]==antibody:
            sam_contr.extend(expand(["{sample}-{control}"], sample = sample, control = samples.loc[sample]["control"]))
    return sam_contr

def get_macs2_peaks_ab(wildcards):
    return expand(
        "results/macs2_callpeak/{sam_contr_peak}_peaks.{peak}Peak",
        sam_contr_peak = get_sample_control_peak_combinations_list_ab(wildcards.antibody), peak =config["params"]["peak-analysis"]
    )

def get_plot_homer_annotatepeaks_input():
    return expand("results/homer/annotate_peaks/{sam_contr_peak}_peaks.annotatePeaks.txt",
        sam_contr_peak = get_sample_control_peak_combinations_list()
    )

def get_samtools_view_filter_input(wildcards):
    return ["results/picard_dedup/{sample}.bam", f"{assembly_path}{assembly}.blacklist.sorted.complement".format(
        prefix="chr{chr}_".format(chr=chromosome) if chromosome else "",
        build=build
    )]

def exists_multiple_groups(antibody):
    return len(samples[samples["antibody"] == antibody]["group"].unique()) > 1

def exists_replicates(antibody):
    return len(samples[samples["antibody"] == antibody]["sample"].unique()) > 1

def not_all_control(antibody):
    return not all([is_control(sample) for sample in samples[samples["antibody"] == antibody]["sample"].unique()])

def get_controls_of_antibody(antibody):
    groups = samples[samples["antibody"] == antibody]["group"]
    controls = samples[pd.isnull(samples["control"])]
    return controls[controls["group"].index.isin(list(groups.index))]["sample"]

def get_samples_of_antibody(antibody):
    groups = samples[samples["antibody"] == antibody]["group"]
    treated = samples[pd.notnull(samples["control"])]
    return treated[treated["group"].index.isin(list(groups.index))]["sample"]

def get_map_reads_input(wildcards):
    if is_sra_pe(wildcards.sample, wildcards.unit):
        return ["results/trimmed/{sample}-{unit}_1.fastq.gz", "results/trimmed/{sample}-{unit}_2.fastq.gz"]
    elif is_single_end(wildcards.sample, wildcards.unit):
        return "results/trimmed/{sample}-{unit}.fastq.gz"
    return ["results/trimmed/{sample}-{unit}_1.fastq.gz", "results/trimmed/{sample}-{unit}_2.fastq.gz"]

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}-{unit}\tSM:{sample}-{unit}\tPL:{platform}'".format(
        sample=wildcards.sample,
        unit=wildcards.unit,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])

def get_is_bigwig(wc):
    big_wigs = []
    for sample in samples.index:
        if not is_control(sample):
            big_wigs.extend(expand(["results/big_wig/{sample}-{control}_subtracted.bigWig"], sample = sample, control = samples.loc[sample]["control"]))
    return big_wigs

def get_multiqc_input(wildcards):
    multiqc_input = []
    for (sample, unit) in units.index:
        reads = [ "1", "2" ]
        if config["params"]["trimming"]["activate"]:
            if is_sra_pe(sample, unit):
                multiqc_input.extend(expand (["logs/cutadapt/{sample}-{unit}.pe.log"],
                sample = sample, unit = unit))
            elif is_single_end(sample, unit):
                reads = [ "0" ]
                multiqc_input.extend(expand (["logs/cutadapt/{sample}-{unit}.se.log"],
                sample = sample, unit = unit))
            else:
                multiqc_input.extend(expand (["logs/cutadapt/{sample}-{unit}.pe.log"],
                sample = sample, unit = unit))
        if not is_sra_pe(sample, unit) and is_single_end(sample,unit):
            reads = [ "0" ]
        multiqc_input.extend(
            expand (
                [
                    "results/qc/fastqc/{sample}.{unit}.{reads}_fastqc.zip",
                    "results/qc/fastqc/{sample}.{unit}.{reads}.html",
                    "results/mapped/{sample}-{unit}.mapped.flagstat",
                    "results/mapped/{sample}-{unit}.mapped.idxstats",
                    "results/mapped/{sample}-{unit}.mapped.stats.txt"
                ],
                sample = sample,
                unit = unit,
                reads = reads
            )
        )
        if config["params"]["trimming"]["activate"]:
            multiqc_input.extend(
                expand (
                    [
                        "results/qc/fastqc/trimmed_{sample}.{unit}.{reads}_fastqc.zip",
                        "results/qc/fastqc/trimmed_{sample}.{unit}.{reads}.html"
                    ],
                    sample = sample,
                    unit = unit,
                    reads = reads                  
                )
            )
    for sample in samples.index:
        multiqc_input.extend(
            expand (
                [
                    "results/picard_dedup/{sample}.metrics.txt",
                    "results/picard_dedup/{sample}.picard_dedup.flagstat",
                    "results/picard_dedup/{sample}.picard_dedup.idxstats",
                    "results/picard_dedup/{sample}.picard_dedup.stats.txt",
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.flagstat",
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.idxstats",
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.stats.txt",
                    "results/phantompeakqualtools/{sample}.phantompeak.spp.out",
                    "results/phantompeakqualtools/{sample}.spp_correlation_mqc.tsv",
                    "results/phantompeakqualtools/{sample}.spp_nsc_mqc.tsv",
                    "results/phantompeakqualtools/{sample}.spp_rsc_mqc.tsv"
                ],
                sample = sample
            )
        )
        if config["params"]["deeptools-plots"]["activate"]:
            multiqc_input.extend(
                expand(
                    [
                        "results/deeptools/plot_profile_data_{peak}.tab"
                    ],
                    sample=sample,
                    peak = config["params"]["peak-analysis"]
                )
            )

        if not is_control(sample):
            multiqc_input.extend(
                expand (
                    [
                        "results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt",
                        "results/deeptools/{sample}-{control}.fingerprint_counts.txt",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.xls"
                    ],
                sample = sample,
                control = samples.loc[sample]["control"],
                peak = config["params"]["peak-analysis"]
            )
        )
        if config["params"]["lc_extrap"]["activate"]:
                multiqc_input.extend( expand(["results/preseq/{sample}.lc_extrap"], sample = sample))
        if config["params"]["picard_metrics"]["activate"]:
            if not config["single_end"]:
                multiqc_input.extend(
                    expand(
                        [
                            "results/qc/multiple_metrics/{sample}.insert_size_metrics",
                            "results/qc/multiple_metrics/{sample}.insert_size_histogram.pdf",
                        ],
                    sample = sample
                )
            )
            multiqc_input.extend(
                expand (
                    [                        
                        "results/qc/multiple_metrics/{sample}.alignment_summary_metrics",
                        "results/qc/multiple_metrics/{sample}.base_distribution_by_cycle_metrics",
                        "results/qc/multiple_metrics/{sample}.base_distribution_by_cycle.pdf",
                        "results/qc/multiple_metrics/{sample}.quality_by_cycle_metrics",
                        "results/qc/multiple_metrics/{sample}.quality_by_cycle.pdf",
                        "results/qc/multiple_metrics/{sample}.quality_distribution_metrics",
                        "results/qc/multiple_metrics/{sample}.quality_distribution.pdf"
                    ], 
                sample = sample
            )
        )
    return multiqc_input

def all_input(wildcards):
    do_annot = config["params"]["peak-annotation-analysis"]["activate"]
    do_peak_qc = config["params"]["peak-qc"]["activate"]
    do_consensus_peak = config["params"]["consensus-peak-analysis"]["activate"]

    wanted_input = []

    # QC with fastQC and multiQC
    wanted_input.extend([
        "results/qc/multiqc/multiqc.html"
    ])

    # trimming reads
    if config["params"]["trimming"]["activate"]:
        for (sample, unit) in units.index:
            if is_sra_pe(sample, unit):
                wanted_input.extend(
                    expand (
                        [
                            "results/trimmed/{sample}-{unit}_1.fastq.gz",
                            "results/trimmed/{sample}-{unit}_2.fastq.gz",
                            "results/trimmed/{sample}-{unit}.pe.qc.txt"
                        ],
                        sample = sample,
                        unit = unit
                )
            )
            elif is_single_end(sample, unit):
                wanted_input.extend(expand(
                        [
                            "results/trimmed/{sample}-{unit}.fastq.gz",
                            "results/trimmed/{sample}-{unit}.se.qc.txt"
                        ],
                        sample = sample,
                        unit = unit
                    )
                )
            else:
                wanted_input.extend(
                    expand (
                        [
                            "results/trimmed/{sample}-{unit}_1.fastq.gz",
                            "results/trimmed/{sample}-{unit}_2.fastq.gz",
                            "results/trimmed/{sample}-{unit}.pe.qc.txt"
                        ],
                        sample = sample,
                        unit = unit
                )
            )

    # mapping, merging and filtering bam-files
    for sample in samples.index:
        wanted_input.extend(
            expand (
                [
                    "results/phantompeakqualtools/{sample}.phantompeak.pdf"
                ],
                sample = sample
            )
        )
        
        if config["params"]["lc_extrap"]["activate"]:
                wanted_input.extend( expand(["results/preseq/{sample}.complexity_measures"], sample = sample))

        if config["params"]["deeptools-plots"]["activate"]:
            wanted_input.extend(
                expand(
                    [
                        "results/deeptools/plot_profile_{peak}.pdf",
                        "results/deeptools/heatmap_{peak}.pdf",
                        "results/deeptools/heatmap_matrix_{peak}.tab"
                    ],
                    sample=sample,
                    peak = config["params"]["peak-analysis"]
                )
            )

        if not is_control(sample):
            with checkpoints.get_gsize.get().output[0].open() as f:
                # only produce the following files, if a gsize is specified
                if f.read().strip() != "":
                    if do_annot:
                        wanted_input.extend(
                            expand(
                                [
                                    "results/homer/annotate_peaks/{sample}-{control}.{peak}_peaks.annotatePeaks.txt"
                                ],
                                sample = sample,
                                control = samples.loc[sample]["control"],
                                peak = config["params"]["peak-analysis"]
                            )
                        )
                        if do_peak_qc:
                            wanted_input.extend(
                                expand(
                                    [
                                        "results/homer/plots/plot_{peak}_annotatepeaks_summary.txt",
                                        "results/homer/plots/plot_{peak}_annotatepeaks_summary.pdf"
                                    ],
                                    peak = config["params"]["peak-analysis"]
                                )
                            )
                    if do_consensus_peak:
                        for antibody in samples["antibody"]:
                            if (exists_multiple_groups(antibody) or exists_replicates(antibody)) and not_all_control(antibody):
                                wanted_input.extend(
                                    expand(
                                        [
                                            "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.saf",
                                            "results/macs2_merged_expand/plots/{antibody}.consensus_{peak}-peaks.boolean.intersect.plot.pdf",
                                            "results/IGV/consensus/merged_library.{antibody}.consensus_{peak}-peaks.igv.txt"
                                        ],
                                        peak = config["params"]["peak-analysis"],
                                        antibody = antibody
                                    )
                                )
                                if do_annot:
                                    wanted_input.extend(
                                        expand(
                                            [
                                                "results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.txt",
                                                "results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.boolean.annotatePeaks.txt",
                                                "results/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts",
                                                "results/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts.summary",
                                                "results/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts.jcounts",
                                                "results/deseq2/dss_rld/{antibody}.consensus_{peak}-peaks.dds.rld.RData",
                                                "results/deseq2/plots/{antibody}.consensus_{peak}-peaks.pca_plot.pdf",
                                                "results/deseq2/plots/{antibody}.consensus_{peak}-peaks.heatmap_plot.pdf",
                                                "results/deseq2/pca_vals/{antibody}.consensus_{peak}-peaks.pca.vals.txt",
                                                "results/deseq2/dists/{antibody}.consensus_{peak}-peaks.sample.dists.txt",
                                                "results/deseq2/sizeFactors/{antibody}.consensus_{peak}-peaks.sizeFactors.RData",
                                                "results/deseq2/sizeFactors/{antibody}.consensus_{peak}-peaks.sizeFactors.sizeFactor.txt",
                                                "results/deseq2/results/{antibody}.consensus_{peak}-peaks.deseq2_results.txt",
                                                "results/deseq2/FDR/{antibody}.consensus_{peak}-peaks.deseq2.FDR_0.01.results.txt",
                                                "results/deseq2/FDR/{antibody}.consensus_{peak}-peaks.deseq2.FDR_0.05.results.txt",
                                                "results/deseq2/FDR/{antibody}.consensus_{peak}-peaks.deseq2.FDR_0.01.results.bed",
                                                "results/deseq2/FDR/{antibody}.consensus_{peak}-peaks.deseq2.FDR_0.05.results.bed",
                                                "results/deseq2/plots/FDR/{antibody}.consensus_{peak}-peaks_FDR_0.01_MA_plot.pdf",
                                                "results/deseq2/plots/FDR/{antibody}.consensus_{peak}-peaks_FDR_0.05_MA_plot.pdf",
                                                "results/deseq2/plots/FDR/{antibody}.consensus_{peak}-peaks_FDR_0.01_volcano_plot.pdf",
                                                "results/deseq2/plots/FDR/{antibody}.consensus_{peak}-peaks_FDR_0.05_volcano_plot.pdf",
                                                "results/deseq2/plots/{antibody}.consensus_{peak}-peaks_sample_corr_heatmap.pdf",
                                                "results/deseq2/plots/{antibody}.consensus_{peak}-peaks_scatter_plots.pdf"
                                            ],
                                            peak = config["params"]["peak-analysis"],
                                            antibody = antibody
                                        )
                                    )
            wanted_input.extend(
                expand(
                    [
                        "results/deeptools/{sample}-{control}.plot_fingerprint.pdf",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_treat_pileup.bdg",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_control_lambda.bdg",
                        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
                        "results/IGV/macs2_callpeak-{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_count.pdf",
                        "results/big_wig/{sample}-{control}_subtracted.bigWig",
                    ],
                    sample = sample,
                    control = samples.loc[sample]["control"],
                    peak = config["params"]["peak-analysis"],
                    antibody = samples.loc[sample]["antibody"]
                )
            )
            if do_peak_qc:
                wanted_input.extend(expand(
                    [
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_frip_score.pdf",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_macs2.pdf"
                    ],
                    peak = config["params"]["peak-analysis"]
                )
            )
            if config["params"]["peak-analysis"] == "broad":
                wanted_input.extend(
                    expand(
                        [
                            "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.gappedPeak"
                        ],
                        sample = sample,
                        control = samples.loc[sample]["control"],
                        peak = config["params"]["peak-analysis"]
                    )
                )
            if config["params"]["peak-analysis"] == "narrow":
                wanted_input.extend(
                    expand(
                        [
                            "results/macs2_callpeak/{sample}-{control}.{peak}_summits.bed"
                        ],
                        sample = sample,
                        control = samples.loc[sample]["control"],
                        peak = config["params"]["peak-analysis"]
                    )
                )
    return wanted_input


genecode = {
    "GRCh38" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz",
        "rmsk" : "hg38"
    },
    "GRCh37" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz",
        "rmsk" : "hg19"
    },
    "GRCm39" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.primary_assembly.annotation.gtf.gz",
        "rmsk" : "mm39"
    },
    "GRCm38" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz",
        "rmsk" : "mm10"
    }
}

genecode_assembly = False
if assembly in genecode.keys():
    genecode_assembly = True

def get_assembly_rmsk(wc):
    if genecode_assembly:
        return genecode[assembly]["rmsk"]
    return assembly