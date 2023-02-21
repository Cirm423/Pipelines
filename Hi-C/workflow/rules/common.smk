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

##Add naming helper functions

#Converts kb and mb notation to base pairs
def bases_to_bp(bases):
    bases = bases.lower()
    if type(bases) == int:
        return str(bases)
    if "mb" in bases:
        return str(int(bases.split("mb")[0]) * 1000000)
    elif "kb" in bases:
        return str(int(bases.split("kb")[0]) * 1000)
    else:
        return bases

#Extract chromosome, start and end from coordinates in config
def split_coords(coords):
    if ":" in coords:
        chrom, coord = coords.split(":")
        start, end = coord.split("-")
        return [chrom,bases_to_bp(start),bases_to_bp(end)]
    else:
        return [coords]

#Make the name to use for enzymes and fragments files

enzyme_file = "_".join(config["params"]["fanc"]["enzyme"].split(","))
fragments_file = "_".join(config["params"]["fanc"]["chr"].split(",")) if config["params"]["fanc"]["chr"] else "all"
bin_sizes = config["params"]["fanc"]["hic"]["bin_size"].split(",")
analysis_resolution = config["params"]["fanc"]["analysis"]["bin_size"]
TAD_res = bases_to_bp(analysis_resolution)

# Get list of regions for the matrix analysis
if config["params"]["fanc"]["regions"]:
    regions = config["params"]["fanc"]["regions"].split(",")
else:
    regions = False

# #List of groups to remove || No use anymore since the way controls are handled changed
# no_group = ['control','input']
# #We remove control and input group from the groups that will be used to call peaks.
# groups = [x for x in groups if x not in no_group]
##### wildcard constraints #####

wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"]),
    group = "|".join(groups),
    sample_group = "|".join(samples.index) + "|" + "|".join(groups)

#Check that settings that allow strings have valid values
if config["params"]["fanc"]["filter"]["multimap"]:
    assert config["params"]["fanc"]["filter"]["multimap"]==True or config["params"]["fanc"]["filter"]["multimap"] == "strict", "multimap setting must be True, False or strict"

#Used only for fanc directionality
# if config["params"]["fanc"]["analysis"]["activate"] and not config["params"]["fanc"]["analysis"]["pca_only"]:
#     tad_form = config["params"]["fanc"]["analysis"]["TAD_format"]
#     assert tad_form == "bed" or tad_form == "gff" or tad_form == "bigwig", "The TAD output format must be either bed gff or bigwig, please check this option"

#Check that analysis bin size is in matrix bin sizes if analysis is activated
if config["params"]["fanc"]["analysis"]["activate"]:
    assert analysis_resolution in bin_sizes, "The analysis bin size must be one of the bin sizes specified in the hic section"

##### reference genomes #####

genecode = {
    "GRCh38-PRI" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz",
        "rmsk" : "hg38"
    },
    "GRCh38-CHR" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz",
        "rmsk" : "hg38"
    },
    "GRCh37" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/gencode.v43lift37.annotation.gtf.gz",
        "rmsk" : "hg19"
    },
    "GRCm39-PRI" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz",
        "rmsk" : "mm39"
    },
    "GRCm39-CHR" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz",
        "rmsk" : "mm39"
    },
    "GRCm38-PRI" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz",
        "rmsk" : "mm10"
    },
    "GRCm38-CHR" : {
        "assembly" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz",
        "gtf" : "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
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

def get_tadlib_input(wildcards):
    sample_g = samples[samples['group'] == wildcards.group]
    return expand(["results/cooler/{sample}.{enzyme}.{fragments}-{resolution}.mcool"],
        sample = sample_g["sample"].index,
        enzyme = enzyme_file,
        fragments = fragments_file,
        resolution = analysis_resolution
    )

def get_pairs_files(wildcards):
    if config["params"]["fanc"]["merge_groups"]:
        sample_g = samples[samples['group'] == wildcards.sample_group]
        return expand([f"results/pairs/{{sample}}.{enzyme_file}.{fragments_file}.pairs"],
            sample = sample_g["sample"].index
        )
    else:
        return f"results/pairs/{{sample_group}}.{enzyme_file}.{fragments_file}.pairs"

def get_hic_files(wildcards):
    if config["params"]["fanc"]["merge_groups"]:
        sample_g = samples[samples['group'] == wildcards.sample_group]
        return expand([f"results/hic/{{sample_group}}.{enzyme_file}.{fragments_file}.hic"],
            sample_group = sample_g["sample"].index
        )
    else:
        return f"results/pairs/{{sample_group}}.{enzyme_file}.{fragments_file}.pairs"

def get_coord_params(wildcards):
    coords = split_coords(wildcards.region)
    if len(coords) == 1:
        return f"-C {coords[0]}"
    elif len(coords) == 3:
        return f"-C {coords[0]} -S {coords[1]} -E {coords[2]}"

#Get the group of the sample and return the hitad file for the group
def get_hitad_output(wildcards):
    group_s = samples.loc[wildcards.sample]["group"]
    return expand(["results/hitad/{group}.{enzyme}.{fragments}-{resolution}.hitad.txt"],
        group = group_s,
        enzyme = enzyme_file,
        fragments = fragments_file,
        resolution = analysis_resolution
    )

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
                    "results/mapped/{sample}_R{reads}.flagstat",
                    "results/mapped/{sample}_R{reads}.idxstats",
                    "results/mapped/{sample}_R{reads}.stats.txt"
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

    return multiqc_input

def all_input(wildcards):
    merge_groups = config["params"]["fanc"]["merge_groups"]
    juicer = config["params"]["fanc"]["hic"]["juicer"]
    do_analysis = config["params"]["fanc"]["analysis"]["activate"]
    multi_TADs = config["params"]["fanc"]["analysis"]["hierarchical_TADs"]

    wanted_input = []

    # QC with fastQC and multiQC, individual stuff
    wanted_input.extend(
        expand([
        "results/qc/multiqc/multiqc.html",
        "resources/{assembly}.{enzyme}.{fragments}.fragments.bed"
        ],
        assembly = assembly,
        enzyme = enzyme_file,
        fragments = fragments_file
        ))

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
                "results/pairs/{sample}.{enzyme}.{fragments}.pairs"
            ],
            sample = sample,
            enzyme = enzyme_file,
            fragments = fragments_file
        ))

        #Output at sample level
        if not merge_groups:
            wanted_input.extend(expand(
                [
                    "results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic",
                    "results/cooler/{sample_group}.{enzyme}.{fragments}-{resolution}.mcool"
                ],
                sample_group = sample,
                enzyme = enzyme_file,
                fragments = fragments_file,
                resolution = bin_sizes
            ))

            if juicer:
                wanted_input.extend(expand(
                        [
                            "results/juicer/{sample_group}.{enzyme}.{fragments}.juicer.hic",
                        ],
                        sample_group = sample,
                        enzyme = enzyme_file,
                        fragments = fragments_file
                    )
                )

            if do_analysis:
                wanted_input.extend(expand(
                    [
                        "results/pca/matrix.{enzyme}.{fragments}-{resolution}.pca_plot.pdf",
                        "results/matrix_analysis/{sample_group}_{chr}.{enzyme}.{fragments}-{resolution}.distance_decay.pdf"
                    ],
                    sample_group = sample,
                    chr = config["params"]["fanc"]["analysis"]["expected_params"],
                    enzyme = enzyme_file,
                    fragments = fragments_file,
                    resolution = analysis_resolution
                ))
        else:
            #Account for TADlib output that needs single sample files
            if do_analysis and multi_TADs:
                wanted_input.extend(expand(
                    [
                        "results/hitad/{sample}.{enzyme}.{fragments}-{resolution}.hitad.DIs.BedGraph"
                    ],
                    sample = sample,
                    chr = config["params"]["fanc"]["analysis"]["expected_params"],
                    enzyme = enzyme_file,
                    fragments = fragments_file,
                    resolution = analysis_resolution
                ))
                if regions:
                    for region in regions:
                        #HiTAD plots only work for specific region coordinates, so only use regions that have coords.
                        if len(split_coords(region)) > 1:
                            wanted_input.extend(expand(
                                [
                                    "results/hitad/{sample}.{enzyme}.{fragments}.{region}-{resolution}.hitad.png"
                                ],
                                sample = sample,
                                region = region,
                                enzyme = enzyme_file,
                                fragments = fragments_file,
                                resolution = analysis_resolution
                            ))

    if merge_groups:
        for group in groups:
            wanted_input.extend(expand(
                    [
                        "results/hic/{sample_group}.{enzyme}.{fragments}-{resolution}.hic",
                        "results/cooler/{sample_group}.{enzyme}.{fragments}-{resolution}.mcool"
                    ],
                    sample_group = group,
                    enzyme = enzyme_file,
                    fragments = fragments_file,
                    resolution = bin_sizes
                )
            )

            if juicer:
                wanted_input.extend(expand(
                        [
                            "results/juicer/{sample_group}.{enzyme}.{fragments}.juicer.hic",
                        ],
                        sample_group = group,
                        enzyme = enzyme_file,
                        fragments = fragments_file
                    )
                )

            if do_analysis:
                wanted_input.extend(expand(
                    [
                        "results/matrix_analysis/{group}_{chr}.{enzyme}.{fragments}-{resolution}.distance_decay.pdf",
                        "results/matrix_analysis/loops/{group}.{enzyme}.{fragments}-{resolution}.merged.bedpe",
                    ],
                    group = group,
                    chr = config["params"]["fanc"]["analysis"]["expected_params"],
                    enzyme = enzyme_file,
                    fragments = fragments_file,
                    resolution = analysis_resolution
                ))

                if not multi_TADs:
                    wanted_input.extend(expand(
                        [
                            "results/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.tad.bed",
                            "results/domaincaller/{group}.{enzyme}.{fragments}-{resolution}.DIs.bedgraph"
                        ],
                        group = group,
                        chr = config["params"]["fanc"]["analysis"]["expected_params"],
                        enzyme = enzyme_file,
                        fragments = fragments_file,
                        resolution = analysis_resolution
                    ))
                else:
                    wanted_input.extend(expand(
                        [
                            "results/hitad/{group}.{enzyme}.{fragments}-{resolution}.hitad.txt",
                        ],
                        group = group,
                        chr = config["params"]["fanc"]["analysis"]["expected_params"],
                        enzyme = enzyme_file,
                        fragments = fragments_file,
                        resolution = analysis_resolution
                    ))

                if regions:
                    for region in regions:
                        wanted_input.extend(expand(
                            [
                                "results/matrix_analysis/compartments/{group}.{enzyme}.{fragments}.{region}-{resolution}.compartments.png",
                            ],
                            group = group,
                            region = region,
                            enzyme = enzyme_file,
                            fragments = fragments_file,
                            resolution = analysis_resolution
                        ))

                        if not multi_TADs:
                            #HiTAD plots only work for specific region coordinates, so only use regions that have coords.
                            if len(split_coords(region)) > 1:
                                wanted_input.extend(expand(
                                    [
                                        "results/domaincaller/{group}.{enzyme}.{fragments}.{region}-{resolution}.tad.png"
                                    ],
                                    group = group,
                                    region = region,
                                    enzyme = enzyme_file,
                                    fragments = fragments_file,
                                    resolution = analysis_resolution
                                ))

    return wanted_input