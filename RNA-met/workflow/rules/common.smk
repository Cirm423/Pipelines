from snakemake.utils import validate
import pandas as pd
import os
from smart_open import open
import yaml
import re

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

#Get the proper star folder name
if config["single_end"]:
    star_lib = "se"
else:
    star_lib = "pe"

#Check that the control samples are actually in the sample table

# samples actually used as controls by other samples
control_samples = set()
for s in samples.index:
    ctrl = samples.loc[s, "control"]
    if pd.notna(ctrl) and str(ctrl).strip():
        control_samples.add(str(ctrl).strip())

# validate
assert control_samples.issubset(set(samples.index)), \
    f"Unknown controls: {control_samples - set(samples.index)}"

bad = [c for c in control_samples if pd.notna(samples.loc[c, "control"]) and str(samples.loc[c, "control"]).strip()]
if bad:
    raise ValueError(f"Samples used as controls cannot themselves have a control value: {bad}")


##### wildcard constraints #####

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in samples.index),
    control = "|".join(re.escape(c) for c in (set(samples["control"].dropna().astype(str).str.strip()) | {"no_control"}))

####### helpers ###########

# treatment samples are everything that is not referenced as a control
treatment_samples = [s for s in samples.index if s not in control_samples]

def sample_control(sample):
    """Return the actual control sample name, or 'no_control'."""
    ctrl = samples.loc[sample, "control"]
    return str(ctrl).strip() if pd.notna(ctrl) and str(ctrl).strip() else "no_control"

def peak_label(sample):
    ctrl = sample_control(sample)
    return f"{sample}-{ctrl}"

def peak_label_with_type(sample, peak_type):
    return f"{peak_label(sample)}.{peak_type}"

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

def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra_accession"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    if ((all_paired & ~sra_null)).all():
        if config["single_end"]:
            all_paired = (~paired).all()
    return all_paired


def has_only_sra_accession(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq1"]) and pd.isnull(units.loc[(sample, unit), "fq2"]) \
           and not pd.isnull(units.loc[(sample, unit), "sra_accession"])

def is_sra_se(sample, unit):
    return has_only_sra_accession(sample, unit) and config["single_end"]

def is_sra_pe(sample, unit):
    return has_only_sra_accession(sample, unit) and not config["single_end"]

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra_accession"]
        if is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("sra-pe-reads/{accession}_{read}.fixed.fastq.gz", accession=accession, read=[1, 2])
        if is_sra_se(wildcards.sample, wildcards.unit):
            return f"sra-se-reads/{accession}.fixed.fastq.gz"

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )

def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files

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

def get_fastqs(wc):
    unit = units.loc[wc.sample]
    if wc.read != "single":
        fq = "fq{}".format(wc.read[-1])
    else:
        fq = "fq1"
    if not units.loc[wc.sample, fq][0].endswith("gz"):
        return units.loc[wc.sample, fq].tolist()

def get_fastqs_gz(wc):
    if config["params"]["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}_{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    unit = units.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra_accession"]
        if not config["single_end"]:
            return expand(
                "sra-pe-reads/{accession}_{read}.fastq.gz", accession=accession, read=wc.read[-1]
            )
        if config["single_end"]:
            return expand(
                "sra-se-reads/{accession}.fastq.gz", accession=accession
            )
    if wc.read != "single":
        fq = "fq{}".format(wc.read[-1])
    else:
        fq = "fq1"
    out_list = units.loc[wc.sample, fq].tolist()
    if units.loc[wc.sample, fq][0].endswith("gz"):
        return out_list


def get_map_reads_input_R1(wildcards):
    unit=units.loc[wildcards.sample]
    if all(pd.isna(unit["fq1"])):
        if not config["single_end"]:
            return "results/merged/{sample}_R1.fastq.gz"
        if config["single_end"]:
            return "results/merged/{sample}_single.fastq.gz"
    ext = units.loc[wildcards.sample]["fq1"][0]
    if ext.endswith("gz") or config["params"]['trimming']['activate']:
        ending = ".gz"
    else:
        ending = ""
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R1.fastq" + ending
    return "results/merged/{sample}_single.fastq" + ending


def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        unit=units.loc[wildcards.sample]
        if all(pd.isna(unit["fq1"])):
            if not config["single_end"]:
                return "results/merged/{sample}_R2.fastq.gz"
            if config["single_end"]:
                return ""            
        ext = units.loc[wildcards.sample]["fq1"][0]
        if ext.endswith("gz") or config["params"]['trimming']['activate']:
            ending = ".gz"
        else:
            ending = ""
        return ("results/merged/{sample}_R2.fastq" + ending,)
    return ""

def macs2_control_input(wildcards):
    ctrl = samples.loc[wildcards.sample, "control"]
    if pd.isna(ctrl) or str(ctrl).strip() == "":
        return []
    return f"results/filtered/{star_lib}/{wildcards.control}.filtered.sortedByCoord.out.bam"

#Old that assumes there is an input
# def is_control(sample):
#     control = samples.loc[sample]["control"]
#     return pd.isna(control) or pd.isnull(control)

def get_sample_control_peak_combinations_list():
    return [
        f"{sample}-{sample_control(sample)}.{config["params"]["peak-analysis"]}"
        for sample in treatment_samples
    ]
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
    return [
        peak_label_with_type(s, config["params"]["peak-analysis"])
        for s in treatment_samples
        if samples.loc[s, "antibody"] == antibody
    ]

def get_sample_control_peak_combinations_list_ab_nop(antibody):
    return [
        peak_label(s)
        for s in treatment_samples
        if samples.loc[s, "antibody"] == antibody
    ]

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

#Changed this to consider inputs having the same antibody as their samples not counting for consensus
def exists_replicates(antibody):
    return sum(
        1 for s in treatment_samples
        if samples.loc[s, "antibody"] == antibody
    ) > 1

def not_all_control(antibody):
    return any(
        s in treatment_samples and samples.loc[s, "antibody"] == antibody
        for s in samples.index
    )

def get_controls_of_antibody(antibody):
    return [
        c
        for c in control_samples
        if samples.loc[c, "antibody"] == antibody
    ]

def get_samples_of_antibody(antibody):
    return [
        s
        for s in treatment_samples
        if samples.loc[s, "antibody"] == antibody
    ]
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
    bigwigs = []
    for s in treatment_samples:
        ctrl = sample_control(s)
        if ctrl == "no_control":
            bigwigs.append(f"results/big_wig/{s}-no_control.bigWig")
        else:
            bigwigs.append(f"results/big_wig/{s}-{ctrl}_subtracted.bigWig")
    return bigwigs
    
def bigwig_path_for(sample):
    label = peak_label(sample)
    ctrl = sample_control(sample)
    if ctrl == "no_control":
        return f"results/big_wig/{label}.bigWig"
    return f"results/big_wig/{label}_subtracted.bigWig"

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
                    "results/qc/fastqc/{sample}.{unit}.{reads}.html"
                ],
                star_lib = star_lib,
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
                    "results/filtered/{star_lib}/{sample}.metrics.txt",
                    "results/filtered/{star_lib}/{sample}.filtered.flagstat",
                    "results/filtered/{star_lib}/{sample}.filtered.idxstats",
                    "results/filtered/{star_lib}/{sample}.filtered.stats.txt",
                ],
                star_lib = star_lib,
                sample = sample
            )
        )
        if config["params"]["phantompeaks"]["activate"]:
            multiqc_input.extend(
                expand(
                    [
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

        if sample in treatment_samples:
            multiqc_input.extend(
                expand (
                    [
                        "results/macs2_callpeak/{label}.{peak}_peaks.xls"
                    ],
                sample = sample,
                label = peak_label(sample),
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
    do_deseq2 = config["params"]["consensus-peak-analysis"]["deseq2"]["activate"]

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
        if config["params"]["phantompeaks"]["activate"]:
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

        if sample in treatment_samples:
            label = peak_label(sample)
            if do_annot:
                wanted_input.extend(
                    expand(
                        [
                            "results/homer/annotate_peaks/{label}.{peak}_peaks.annotatePeaks.txt",
                            "results/homer/motifs/{label}.{peak}_peaks/homerMotifs.all.motifs"
                        ],
                        sample = sample,
                        label = label,
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
                    if exists_replicates(antibody):
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
                                        "results/homer/motifs_consensus/{antibody}.consensus_{peak}-peaks/homerMotifs.all.motifs",                                              
                                        ],
                                        peak = config["params"]["peak-analysis"],
                                        antibody = antibody
                                    )
                                )
                        if do_deseq2:
                            if (exists_multiple_groups(antibody) or exists_replicates(antibody)) and not_all_control(antibody):
                                wanted_input.extend(
                                    expand(
                                        [                        
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
                        "results/macs2_callpeak/{label}.{peak}_treat_pileup.bdg",
                        "results/macs2_callpeak/{label}.{peak}_control_lambda.bdg",
                        "results/macs2_callpeak/{label}.{peak}_peaks.{peak}Peak",
                        "results/IGV/macs2_callpeak-{peak}/merged_library.{label}.{peak}_peaks.igv.txt",
                        "results/macs2_callpeak/plots/plot_{peak}_peaks_count.pdf",
                        "{bigwig}",
                    ],
                    label = label,
                    peak = config["params"]["peak-analysis"],
                    bigwig = bigwig_path_for(sample)
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
                            "results/macs2_callpeak/{label}.{peak}_peaks.gappedPeak"
                        ],
                        label = peak_label(sample),
                        peak = config["params"]["peak-analysis"]
                    )
                )
            if config["params"]["peak-analysis"] == "narrow":
                wanted_input.extend(
                    expand(
                        [
                            "results/macs2_callpeak/{label}.{peak}_summits.bed"
                        ],
                        label = peak_label(sample),
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