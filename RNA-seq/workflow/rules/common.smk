import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

if config["params"]["pca"]["activate"]:
    assert(config["params"]["pca"]["activate"] and config["params"]["diffexp"]["activate"]), "Pca cannot be activated without activating diffexp"
assert(not (config["params"]["diffexp"]["activate"] and (config["params"]["single"]["activate"] or config["params"]["TE_single"]["activate"]))), "Single or TE_single cannot be activated at the same time as diffexp"

alternative_types = ["greater","less","two.sided"]
assert  config["params"]["single"]["alternative"] in alternative_types and config["params"]["TE_single"]["alternative"] in alternative_types, "The single mode alternative must be either greater, less or two.sided"

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str, "condition": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

assembly = config['resources']["ref"]["assembly"]
assembly_path = config['resources']['path'] + config['resources']['ref']['assembly'] + "/"

if config["params"]["single"]["activate"] or config["params"]["TE_single"]["activate"]:
    assert(len(samples.loc[samples["condition"]=="control"])==1), "For single analysis one of the samples has to be called 'control'"

wildcard_constraints:
    sample = "|".join(samples.index)

#Get the proper star folder name
if config["single_end"]:
    lib_end = "se"
else:
    lib_end = "pe"

if config["params"]["2nd_pass"]["activate"]:
    num = "2"
else:
    num = ""

star_lib = lib_end + num

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
if config['resources']['ref']['assembly'] in genecode.keys():
    genecode_assembly = True

def get_assembly_rmsk(wc):
    if genecode_assembly:
        return genecode[config['resources']['ref']['assembly']]["rmsk"]
    return config['resources']['ref']['assembly']

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
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


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
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

# Imported from Chipseq

def has_only_sra_accession(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq1"]) and pd.isnull(units.loc[(sample, unit), "fq2"]) \
           and not pd.isnull(units.loc[(sample, unit), "sra"])

def is_sra_se(sample, unit):
    return has_only_sra_accession(sample, unit) and config["single_end"]

def is_sra_pe(sample, unit):
    return has_only_sra_accession(sample, unit) and not config["single_end"]

def get_individual_fastq(wildcards):
    if ( wildcards.read == "0" or wildcards.read == "1" ):
        if is_sra_se(wildcards.sample, wildcards.unit):
            return expand("sra-se-reads/{accession}.fastq.gz",
                            accession=units.loc[ (wildcards.sample, wildcards.unit), "sra" ])
        elif is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("sra-pe-reads/{accession}_1.fastq.gz",
                            accession=units.loc[ (wildcards.sample, wildcards.unit), "sra" ])
        else:
            return units.loc[ (wildcards.sample, wildcards.unit), "fq1" ]
    elif wildcards.read == "2":
        if is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("sra-pe-reads/{accession}_2.fastq.gz",
                          accession=units.loc[ (wildcards.sample, wildcards.unit), "sra" ])
        else:
            return units.loc[ (wildcards.sample, wildcards.unit), "fq2" ]

def get_individual_trimmed_fastq(wildcards):
    if wildcards.read == "0":
        return expand("results/trimmed/{sample}_{unit}_single.fastq.gz",
                    sample = wildcards.sample, unit = wildcards.unit)
    elif wildcards.read == "1":
        return expand("results/trimmed/{sample}_{unit}_R1.fastq.gz",
                        sample = wildcards.sample, unit = wildcards.unit)
    elif wildcards.read == "2":
        return expand("results/trimmed/{sample}_{unit}_R2.fastq.gz",
                        sample = wildcards.sample, unit = wildcards.unit)

# Original from RNA-seq

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

def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = [0.5]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

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
        accession = unit["sra"]
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

def get_fastqs(wc):
    unit = units.loc[wc.sample]
    if wc.read != "single":
        fq = "fq{}".format(wc.read[-1])
    else:
        fq = "fq1"
    if not units.loc[wc.sample, fq][0].endswith("gz"):
        return units.loc[wc.sample, fq].tolist()

def get_contrast(wildcards):
    return config["params"]["diffexp"]["contrasts"][wildcards.contrast]

#To be removed, here only for reference
# def path_merged_cond(path):
#     if config["params"]["mergeReads"]["activate"]:
#         path = path.replace('?','{sample}')
#     else:
#         path = path.replace('?','{sample}-{unit}')
#     if config["params"]["2nd_pass"]["activate"]:
#         path = path.replace('¿',"2")
#     else:
#         path = path.replace('¿',"")
#     return path

# def path_merged_cond_mqc(path):
#     if config["params"]["mergeReads"]["activate"]:
#         return path.replace('?','{unit.sample_name}')
#     else:
#         return path.replace('?','{unit.sample_name}-{unit.unit_name}')

def get_deseq2_end(wc):
    for unit in units.itertuples():
        if is_paired_end(unit.sample_name):
            return "FALSE"
        return "TRUE"

def get_multiqc_input(wildcards):
    multiqc_input = []
    for (sample, unit) in units.index:
        reads = [ "1", "2" ]
        if not is_sra_pe(sample, unit):
            if not is_paired_end(sample):
                reads = [ "0" ]
        multiqc_input.extend(
            expand (
                [
                    "results/qc/fastqc/{sample}.{unit}.{reads}_fastqc.zip",
                    "results/qc/fastqc/{sample}.{unit}.{reads}.html",
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
                        "results/qc/fastqc/trimmed_{sample}.{unit}.{reads}.html",
                    ],
                    sample = sample,
                    unit = unit,
                    reads = reads
                )
            )
    multiqc_input.extend(expand(
        "results/star/{lib_end}/{sample}/Aligned.sortedByCoord.out.bam",
        lib_end=lib_end,
        sample=samples.sample_name
    ))    
    if config["params"]["rseqc"]["activate"]:
        multiqc_input.extend(
            expand([
                    "results/qc/rseqc/{sample}.junctionanno.junction.bed",
                    "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
                    "results/qc/rseqc/{sample}.infer_experiment.txt",
                    "results/qc/rseqc/{sample}.stats.txt",
                    "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
                    "results/qc/rseqc/{sample}.readdistribution.txt",
                    #"results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
                    "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
                    "logs/rseqc/rseqc_junction_annotation/{sample}.log"
                ],
                sample=samples.sample_name
            ),
        )
                
    return multiqc_input

def all_input(wildcards):

    wanted_input = []

    # QC with fastQC and multiQC
    wanted_input.extend([
        "results/qc/multiqc_report.html",
        f"{assembly_path}{assembly}.annotation.bed"
    ])
    
    if config["params"]["diffexp"]["activate"]:
        wanted_input.extend( expand(
            "results/diffexp/{contrast}.diffexp.tsv",
            contrast=config["params"]["diffexp"]["contrasts"],
        ))
        if config["params"]["diffexp"]["TE"]["activate"]:
            wanted_input.extend(expand(
                "results/diffexp/{contrast}.diffexp.TE.tsv",
                contrast=config["params"]["diffexp"]["contrasts"],
            ))
        if config["params"]["pca"]["activate"]:
            wanted_input.append("results/pca.svg")
            if config["params"]["diffexp"]["TE"]["activate"]:
                wanted_input.append("results/TE_pca.svg")
    else:
        wanted_input.extend(expand(
            "results/rsem/{star_lib}/{sample}/mapped.genes.results",
            sample=units["sample_name"],star_lib=star_lib
        ))
        wanted_input.extend(expand(
            "results/rsem/{star_lib}/{sample}/mapped.isoforms.results",
            sample=units["sample_name"],star_lib=star_lib
        ))
    wanted_input.extend(expand("results/browser/{sample}.str1.bw",sample=units["sample_name"]))
    if get_strandedness(units)[0] != 0.5:
        wanted_input.extend(expand("results/browser/{sample}.str2.bw",sample=units["sample_name"]))
    if config["params"]["single"]["activate"]:
        wanted_input.extend(expand("results/single/{sample}_vs_{control}.tsv",
        sample=samples.loc[samples["condition"]!="control"]["sample_name"],
        control=samples.loc[samples["condition"]=="control"]["sample_name"]))
    if config["params"]["TE_single"]["activate"]:
        wanted_input.extend(expand("results/TE_single/{sample}_vs_{control}.tsv",
        sample=samples.loc[samples["condition"]!="control"]["sample_name"],
        control=samples.loc[samples["condition"]=="control"]["sample_name"]))

    return wanted_input
