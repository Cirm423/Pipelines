from snakemake.utils import validate
import pandas as pd
import os
from smart_open import open
import yaml

#FOR TESTING ONLY
# with open("config/config.yaml",'r') as f:
#      config = yaml.safe_load(f)

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
groups = samples["group"].unique()

assembly = config["resources"]["ref"]["assembly"]
assembly_path = config['resources']['path'] + config['resources']['ref']['assembly'] + "/"

#The pipe should always be run in paired end mode

config["single_end"] = False

#Make the name to use for enzymes and fragments files

enzyme_file = "_".join(config["params"]["fanc"]["enzyme"].split(","))
fragments_file = "_".join(config["params"]["fanc"]["chr"].split(",")) if config["params"]["fanc"]["chr"] else "all"

# #List of groups to remove || No use anymore since the way controls are handled changed
# no_group = ['control','input']
# #We remove control and input group from the groups that will be used to call peaks.
# groups = [x for x in groups if x not in no_group]
##### wildcard constraints #####

wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"]),
    group = "|".join(groups)

#Check that settings that allow strings have valid values
if config["params"]["fanc"]["filter"]["multimap"]:
    assert config["params"]["fanc"]["filter"]["multimap"]==True or config["params"]["fanc"]["filter"]["multimap"] == "strict", "multimap setting must be True, False or strict"

##### reference genomes #####

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

#Not used in ATAC since there are not different antybodies, but may be useful for something

# def is_control(sample):
#     control = samples.loc[sample]["control"]
#     return pd.isna(control) or pd.isnull(control)

# def get_sample_control_peak_combinations_list():
#     sam_contr = []
#     for sample in samples.index:
#         if not is_control(sample):
#             sam_contr.extend(expand(["{sample}-{control}.{peak}"], sample = sample, control = samples.loc[sample]["control"], peak = config["params"]["peak-analysis"]))
#     return sam_contr

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

#Change the function bellow to get the ATAC peaks, not macs2
# def get_macs2_peaks():
#     return expand(
#         "results/macs2_callpeak/{sam_contr_peak}_peaks.{peak}Peak",
#         sam_contr_peak = get_sample_control_peak_combinations_list(), peak =config["params"]["peak-analysis"]
#     )

def get_plot_homer_annotatepeaks_input():
    return expand("results/homer/annotate_peaks/{sam_contr_peak}_peaks.annotatePeaks.txt",
        sam_contr_peak = get_sample_control_peak_combinations_list()
    )

# No use for this one, save just in case for now
# def exists_multiple_groups(antibody):
#     return len(samples[samples["antibody"] == antibody]["group"].unique()) > 1

#Not tested
def exists_replicates(group):
    return len(samples[samples["group"] == group]["sample"].unique()) > 1

#No controls in Hi-C
# #Changed this to allow controls to be in a different group as the samples, i.e. in the case of an input
# def get_controls_of_group(group):
#     sample_g = samples[samples['group'] == group]
#     #controls = samples[pd.isnull(samples["control"])]
#     #return controls[controls["group"].index.isin(list(sample_g.index))]["sample"]
#     if pd.isnull(sample_g["control"]).all():
#         return ""
#     else:
#         treated = samples[pd.notnull(samples["control"])]
#         return expand(["results/genrich/{control}.sorted.bam"],
#             control = list(pd.unique(treated[treated["group"].index.isin(list(sample_g.index))]["control"]))
#             )

def get_samples_of_group(group):
    #Accounting for groups with no control and groups with controls
    sample_g = samples[samples['group'] == group]
    if pd.isnull(sample_g["control"]).all():
        return expand(["results/genrich/{sample}.sorted.bam"],
            sample = sample_g["sample"].index
        )
    else:
        treated = samples[pd.notnull(samples["control"])]
        return expand(["results/genrich/{sample}.sorted.bam"],
            sample = treated[treated["group"].index.isin(list(sample_g.index))]["sample"]
        )

def get_pairs_files(wildcards):
    if config["params"]["fanc"]["merge_groups"]:
        sample_g = samples[samples['group'] == wildcards.sample_group]
        return expand(["results/pairs/{sample}.pairs"],
            sample = sample_g["sample"].index
        )
    else:
        return "results/pairs/{sample_group}.pairs"

def get_unit_R1_of_sample(wildcards):
    unit_list = []
    for unit in units.loc[wildcards.sample, "unit"]:
        if config["params"]["trimming"]["activate"]:
            if is_sra_pe(wildcards.sample,unit):
                unit_list.append(f"results/trimmed/{wildcards.sample}-{unit}_1.fastq.gz")
            elif is_single_end(wildcards.sample, unit):
                unit_list.append(f"results/trimmed/{wildcards.sample}-{unit}.fastq.gz")
            else:
                unit_list.append(f"results/trimmed/{wildcards.sample}-{unit}_1.fastq.gz")
        else:
            if is_sra_se(wildcards.sample, unit):
                unit_list.append(f"resources/sra-se-reads/{units.loc[ (wildcards.sample, unit), 'sra_accession' ]}.fastq.gz")
            elif is_sra_pe(wildcards.sample, unit):
                unit_list.append(f"resources/sra-pe-reads/{units.loc[ (wildcards.sample, unit), 'sra_accession' ]}_1.fastq.gz")
            elif is_single_end(wildcards.sample, unit):
                unit_list.append(units.loc[ (wildcards.sample, unit), "fq1" ])
            else:
                u = units.loc[ (wildcards.sample, unit), ["fq1", "fq2"] ].dropna()
                unit_list.append(f"{u.fq1}")
    return unit_list

def get_unit_R2_of_sample(wildcards):
    unit_list = []
    for unit in units.loc[wildcards.sample, "unit"]:
        if config["params"]["trimming"]["activate"]:
            if is_sra_pe(wildcards.sample,unit):
                unit_list.append(f"results/trimmed/{wildcards.sample}-{unit}_2.fastq.gz")
            elif is_single_end(wildcards.sample, unit):
                unit_list.append("")
            else:
                unit_list.append(f"results/trimmed/{wildcards.sample}-{unit}_2.fastq.gz")
        else:
            if is_sra_se(wildcards.sample, unit):
                unit_list.append("")
            elif is_sra_pe(wildcards.sample, unit):
                unit_list.append(f"resources/sra-pe-reads/{units.loc[ (wildcards.sample, unit), 'sra_accession' ]}_2.fastq.gz")
            elif is_single_end(wildcards.sample, unit):
                unit_list.append("")
            else:
                u = units.loc[ (wildcards.sample, unit), ["fq1", "fq2"] ].dropna()
                unit_list.append(f"{u.fq2}")
    return unit_list

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}-{unit}\tSM:{sample}-{unit}\tPL:{platform}'".format(
        sample=wildcards.sample,
        unit=wildcards.unit,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])

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
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.flagstat",
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.idxstats",
                    "results/bamtools_filtered/{sample}.sorted.bamtools_filtered.stats.txt",
                ],
                sample = sample
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

    # QC with fastQC and multiQC, individual stuff
    wanted_input.extend([
        "results/qc/multiqc/multiqc.html",
        "results/genrich/plots/plot_narrow_peaks_frip_score.pdf",
        "results/genrich/plots/plot_narrow_peaks_count.pdf",
        "results/IGV/big_wig/merged_library.bigWig.igv.txt"
    ])

    if do_peak_qc:
        wanted_input.extend(
            [
                "results/genrich/plots/plot_narrow_peaks_genrich.pdf"
            ]
        )
    
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
    for sample in samples.index:
        wanted_input.extend(expand(
            [   
                "results/bedtools_intersect/{sample}.intersected.bed",
                "results/bedtools_intersect/{sample}.narrow.peaks_frip.tsv"
            ],
            sample = sample
        ))

        if config["params"]["lc_extrap"]["activate"]:
                wanted_input.extend( expand(["results/preseq/{sample}.complexity_measures"], sample = sample))
                
        if ATACseqQC_act:
            wanted_input.extend(expand(
                [
                    "results/qc/ATACseqQC/{sample}/fragmentSizeDistribution.pdf",
                    "results/qc/ATACseqQC/{sample}/PTscore.pdf",
                    "results/qc/ATACseqQC/{sample}/NFRscore.pdf",
                    "results/qc/ATACseqQC/{sample}/TSSEscore.pdf",
                    "results/qc/ATACseqQC/{sample}/cumulativePercentage.pdf",
                    "results/qc/ATACseqQC/{sample}/featureAligndHeatmap.pdf",
                    "results/qc/ATACseqQC/{sample}/TSS_profile.pdf",
                    "results/qc/ATACseqQC/{sample}/CTCF_footprint.pdf",
                    "results/qc/ATACseqQC/{sample}/CTCF_Vplot.pdf"
                ],
                sample = sample
            ))
    for group in groups:
        wanted_input.extend(expand(
                [
                    "results/genrich/{group}.narrowPeak",
                    "results/genrich/{group}.bed",
                    "results/IGV/genrich_peaks/merged_library.{group}.narrow_peaks.igv.txt"

                ],
                group = group
            )
        )
        if do_annot:
            wanted_input.extend(expand(
                [
                    "results/homer/annotate_peaks/{group}.narrow_peaks.annotatePeaks.txt"
                ],
                group = group
            )
        )
        
    #Need to add more files as things are made, at least until peak calling for now
    return wanted_input

#Only used for testing
def test_input(wildcards):
    wanted_input = []
    wanted_input.extend([
        f"resources/{assembly}.{enzyme_file}.{fragments_file}.fragments.bed"
    ])
    for sample in samples.index:
        wanted_input.extend(expand(
            [   
                "results/pairs/{sample}.pairs"
            ],
            sample = sample
        ))
        if not config["params"]["fanc"]["merge_groups"]:
            wanted_input.extend(expand(
                [   
                    "results/hic/{sample_group}.hic",
                    "results/juicer/{sample_group}.juicer.hic"
                ],
                sample_group = sample
            ))
    
    for group in groups:
        if config["params"]["fanc"]["merge_groups"]:
            wanted_input.extend(expand(
                [
                    "results/hic/{sample_group}.hic",
                    "results/juicer/{sample_group}.juicer.hic"  
                ],
                sample_group = group
            ))
    print(wanted_input)
    return wanted_input