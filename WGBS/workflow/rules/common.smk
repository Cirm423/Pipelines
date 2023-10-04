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

assembly = config["resources"]["ref"]["assembly"]
assembly_path = config['resources']['path'] + config['resources']['ref']['assembly'] + "/"
phage_path = config['resources']['path'] + "phage_lambda/"

# Check that the mode is correct

assert config["params"]["mode"] == "bwameth" or config["params"]["mode"] == "bismark", "The pipeline mode has to be either 'bwameth' or 'bismark'"

# Check treatment is in the groups when diff_meth is activated

if config["params"]["diff_meth"]["activate"]:
    assert "treatment" in list(samples["group"]), "treatment needs to be a group when performing differential methylation analysis"

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

def get_bismark_bams(wildcards):
    if config["single_end"]:
        return expand("results/bismark_mapped/{sample}_se.bam", sample = samples.index)
    else:
        return expand("results/bismark_mapped/{sample}_pe.bam", sample = samples.index)

def get_bams(wildcards):
    if config["params"]["mode"] == "bwameth":
        return "results/picard_dedup/{sample}.bam"
    else:
        if config["single_end"]:
            return "results/bismark_mapped/{sample}_se.bam"
        else:
            return "results/bismark_mapped/{sample}_pe.bam"

def get_dedup_bam(wildcards):
    if config["params"]["mode"] == "bwameth":
        return "results/picard_dedup/{sample}.bam"
    else:
        if config["single_end"]:
            return "results/bismark_mapped/{sample}_se.deduplicated.bam"
        else:
            return "results/bismark_mapped/{sample}_pe.deduplicated.bam"

def get_sample_splitting_reports(wildcards):
    if config["single_end"]:
        return ["results/bismark/meth/{sample}.deduplicated_splitting_report.txt","results/bismark_mapped/{sample}_SE_report.txt"]
    else:
        return "results/bismark/meth/{sample}_pe.deduplicated_splitting_report.txt"

def get_sample_splitting_reports_ln(wildcards):
    if config["single_end"]:
        return expand(["results/bismark_mapped/{sample}.deduplicated_splitting_report.txt", "results/bismark_mapped/{sample}_se_SE_report.txt"], sample = samples.index)
    else:
        return expand("results/bismark_mapped/{sample}_pe.deduplicated_splitting_report.txt", sample = samples.index)

def get_sample_reports_phage(wildcards):
    if config["single_end"]:
        return expand("results/bismark_mapped/{sample}-phage_SE_report.txt", sample = samples.index)
    else:
        return expand("results/bismark_mapped/{sample}-phage_PE_report.txt", sample = samples.index)

def get_methylkit_input(wildcards):
    if config["params"]["mode"] == "bwameth":
        return expand(
            [
                "results/methyldackel/{sample}_CpG.methylKit"
            ],
            sample = samples.index
        )
    else:
        return expand(
            [
                "results/bismark/meth_cpg/{sample}-{end}.bismark.cov.gz",
            ],
            sample = samples.index,
            end = "se" if config["single_end"] else "pe"
        )

def get_deduplicated_bams(wildcards):
    if config["params"]["mode"] == "bwameth":
        return "results/picard_dedup/{sample}.bam"
    else:
        return expand(
            [
                "results/bismark_mapped/{{sample}}_{end}.deduplicated.bam",
            ],
            end = "se" if config["single_end"] else "pe"
        )

def get_bedgraphs(wildcards):
    if wildcards.bin == "1":
        if config["params"]["mode"] == "bwameth":
            return "results/methyldackel/{sample}.sorted.bedGraph"
        else:
            return "results/bismark/{sample}.sorted.bedGraph"
    else:
        return "results/big_wig/{sample}.binned.sorted.bedGraph"


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
                    "results/qualimap/{sample}_qualimap"
                ],
                sample = sample
            )
        )
        if config["params"]["mode"] == "bwameth":
            multiqc_input.extend(
                expand(
                    [
                        "results/mapped/{sample}.mapped.flagstat",
                        "results/mapped/{sample}.mapped.idxstats",
                        "results/mapped/{sample}.mapped.stats.txt",
                        "results/picard_dedup/{sample}.metrics.txt",
                        "results/picard_dedup/{sample}.picard_dedup.flagstat",
                        "results/picard_dedup/{sample}.picard_dedup.idxstats",
                        "results/picard_dedup/{sample}.picard_dedup.stats.txt",
                        "results/methyldackel/{sample}_methyldackel.txt"                 
                    ],
                    sample=sample
                )
            )
        else:
            multiqc_input.extend(
                expand(
                    [
                        "results/qc/bismark/bismark2summary.html",
                        "results/qc/bismark/bismark2summary.txt"
                    ],
                    sample = sample
                )
            )
            if config["single_end"]:
                multiqc_input.extend(
                    expand(
                        [                            
                            "results/bismark_mapped/{sample}_SE_report.txt",
                            "results/bismark_mapped/{sample}.deduplication_report.txt",
                            "results/bismark/meth/{sample}.deduplicated_splitting_report.txt",
                            "results/bismark/meth/{sample}-se.M-bias.txt",
                            "results/qc/bismark/{sample}_se.bismark2report.html",
                        ],
                        sample=sample
                    )
                )
            else:
                multiqc_input.extend(
                    expand(
                        [
                            "results/bismark_mapped/{sample}_PE_report.txt",
                            "results/bismark_mapped/{sample}_pe.deduplication_report.txt",
                            "results/bismark/meth/{sample}_pe.deduplicated_splitting_report.txt",
                            "results/bismark/meth/{sample}-pe.M-bias.txt",
                            "results/qc/bismark/{sample}_pe.bismark2report.html",
                        ],
                        sample=sample
                    )
                )            
        if config["params"]["lc_extrap"]["activate"]:
                multiqc_input.extend( expand(["results/preseq/{sample}.lc_extrap"], sample = sample))

    return multiqc_input

def all_input(wildcards):

    wanted_input = []

    # QC with fastQC and multiQC
    wanted_input.extend([
        "results/qc/multiqc/multiqc.html"
    ])

    if config["params"]["phage"]:
        wanted_input.extend([
            "results/qc/bisulfite_conversion_rate.csv"
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

    # Methylation, depending on bismark or bwameth
    for sample in samples.index:
        wanted_input.extend(expand("results/big_wig/{sample}_meth-perc_{bin}bp-bins.bigWig",sample=sample,bin=config["params"]["bigwig_bins"]))
        if config["params"]["lc_extrap"]["activate"]:
            wanted_input.extend( expand(["results/preseq/{sample}.lc_extrap"], sample = sample))
        if config["params"]["mode"] == "bwameth":
            wanted_input.extend(
                expand(
                    [
                        "results/methyldackel/{sample}_{meth}{suffix}"
                    ],
                    sample = sample,
                    meth = ["CpG","CHG","CHH"] if config["params"]["methyldackel"]["comprehensive"] else ["CpG"],
                    suffix = [".methylKit",".bedGraph"] if config["params"]["diff_meth"]["activate"] else [".bedGraph"]
                )
            ) 
        else:
            if config["single_end"]:
                wanted_input.extend(
                    expand(
                        [
                            "results/qc/bismark/{sample}-se.M-bias_R1.png",
                            "results/bismark/meth/{sample}-se.M-bias.txt",
                            "results/bismark/meth/{sample}.deduplicated_splitting_report.txt",
                            "results/bismark/meth_cpg/{sample}-se.bismark.cov.gz",
                            "results/bismark/meth_cpg/{sample}-se.bedGraph.gz"
                        ],
                        sample=sample
                    )
                )
            else:
                wanted_input.extend(
                    expand(
                        [
                            "results/qc/bismark/{sample}-pe.M-bias_R1.png",
                            "results/qc/bismark/{sample}-pe.M-bias_R2.png",
                            "results/bismark/meth/{sample}-pe.M-bias.txt",
                            "results/bismark/meth/{sample}_pe.deduplicated_splitting_report.txt",
                            "results/bismark/meth_cpg/{sample}-pe.bismark.cov.gz",
                            "results/bismark/meth_cpg/{sample}-pe.bedGraph.gz"
                        ],
                        sample=sample
                    )
                )
        
        #PMDs
        if config["params"]["PMDs"]["activate"]:
            wanted_input.extend(
                expand(
                    [
                        "results/methpipe/PMDs/{sample}.pmd"
                    ],
                    sample = sample
                )
            )
            
    #Differential methylation analysis
    if config["params"]["diff_meth"]["activate"]:
        wanted_input.extend([
            "results/diff_meth/methDB.RDS",
            "results/diff_meth/plots/Sample_correlation.pdf",
            "results/diff_meth/plots/Sample_clustering.pdf",
            "results/diff_meth/plots/PCA_screen.pdf",
            "results/diff_meth/plots/PCA.pdf",
            "results/diff_meth/CpG_hypermethylated_25p.tsv",
            "results/diff_meth/CpG_hypomethylated_25p.tsv",
            "results/diff_meth/CpG_all_methylated_diff_25p.tsv",
            "results/diff_meth/CpG_methylated_by_chr_25p.tsv",
            "results/diff_meth/CpG_methylated_annotation_25p.tsv",
            "results/big_wig/CpG_all_methylated_diff.bigWig"
        ])
            
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