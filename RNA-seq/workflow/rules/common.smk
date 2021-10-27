import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

if config["pca"]["activate"]:
    assert(config["pca"]["activate"] and config["diffexp"]["activate"]), "Pca cannot be activated without activating diffexp"
assert(not (config["diffexp"]["activate"] and (config["single"]["activate"] or config["TE_single"]["activate"]))), "Single or TE_single cannot be activated at the same time as diffexp"


samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
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

if config["single"]["activate"] or config["TE_single"]["activate"]:
    assert(len(samples.loc[samples["condition"]=="control"])==1), "For single analysis one of the samples has to be called 'control'"
    assert(config["mergeReads"]["activate"]), "For single analysis, mergeReads must be activated"


def get_final_output():
    if config["diffexp"]["activate"]:
        final_output = expand(
            "results/diffexp/{contrast}.diffexp.tsv",
            contrast=config["diffexp"]["contrasts"],
        )
        if config["diffexp"]["TE"]["activate"]:
            final_output.extend(expand(
                "results/diffexp/{contrast}.diffexp.TE.tsv",
                contrast=config["diffexp"]["contrasts"],
            ))
        if config["pca"]["activate"]:
            final_output.append("results/pca.svg")
            if config["diffexp"]["TE"]["activate"]:
                final_output.append("results/TE_pca.svg")
    #elif config["single"]["activate"]:
    else:
        if is_activated("2nd_pass"):
            num = "2"
        else:
            num = ""
        if pd.isna(units["fq2"][0]):
            if is_activated("mergeReads"):
                final_output = expand(
                    "results/rsem/se{num}/{sample}/mapped.genes.results",
                    sample=units["sample_name"],num=num
                )
                final_output.extend(expand(
                    "results/rsem/se{num}/{sample}/mapped.isoforms.results",
                    sample=units["sample_name"],num=num
                ))
            else:
                final_output = expand(
                    "results/rsem/se{num}/{sample}-{unit}/mapped.genes.results",
                    sample=units["sample_name"],unit=units["unit_name"], num=num
                )
                final_output.extend(expand(
                    "results/rsem/se{num}/{sample}-{unit}/mapped.isoforms.results",
                    sample=units["sample_name"],unit=units["unit_name"], num=num
                ))
        else:
            if is_activated("mergeReads"):
                final_output = expand(
                    "results/rsem/pe{num}/{sample}/mapped.genes.results",
                    sample=units["sample_name"], num=num
                )
                final_output.extend(expand(
                    "results/rsem/pe{num}/{sample}/mapped.isoforms.results",
                    sample=units["sample_name"], num=num
                ))
            else:
                final_output = expand(
                    "results/rsem/pe{num}/{sample}-{unit}/mapped.genes.results",
                    sample=units["sample_name"],unit=units["unit_name"], num=num
                )
                final_output.extend(expand(
                    "results/rsem/pe{num}/{sample}-{unit}/mapped.isoforms.results",
                    sample=units["sample_name"],unit=units["unit_name"], num=num
                ))
    if is_activated("mergeReads"):
        final_output.extend(expand("results/browser/{sample}.str1.bw",sample=units["sample_name"]))
        if get_strandedness(units)[0] != 0.5:
            final_output.extend(expand("results/browser/{sample}.str2.bw",sample=units["sample_name"]))
    else:
        final_output.extend(expand("results/browser/{sample}-{unit}.str1.bw",sample=units["sample_name"],unit=units["unit_name"]))
        if get_strandedness(units)[0] != 0.5:
            final_output.extend(expand("results/browser/{sample}-{unit}.str2.bw",sample=units["sample_name"],unit=units["unit_name"]))
    if is_activated("single"):
        final_output.extend(expand("results/single/{sample}_vs_{control}.tsv",
        sample=samples.loc[samples["condition"]!="control"]["sample_name"],
        control=samples.loc[samples["condition"]=="control"]["sample_name"]))
    if is_activated("TE_single"):
        final_output.extend(expand("results/TE_single/{sample}_vs_{control}.tsv",
        sample=samples.loc[samples["condition"]!="control"]["sample_name"],
        control=samples.loc[samples["condition"]=="control"]["sample_name"]))
    return list(set(final_output))

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
if config['ref']['assembly'] in genecode.keys():
    genecode_assembly = True

def get_assembly_rmsk(wc):
    if genecode_assembly:
        return genecode[config['ref']['assembly']]["rmsk"]
    return config['ref']['assembly']

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        if is_sra_pe(wildcards.sample, wildcards.unit):
            return expand("sra-pe-reads/{accession}_{read}.fastq.gz", accession=accession, read=[1, 2])
        if is_sra_se(wildcards.sample, wildcards.unit):
            return f"sra-se-reads/{accession}.fastq.gz"

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
        if config["single_end"]["activate"]:
            all_paired = (~paired).all()
    return all_paired

# Imported from Chipseq

def has_only_sra_accession(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq1"]) and pd.isnull(units.loc[(sample, unit), "fq2"]) \
           and not pd.isnull(units.loc[(sample, unit), "sra"])

def is_sra_se(sample, unit):
    return has_only_sra_accession(sample, unit) and config["single_end"]["activate"]

def is_sra_pe(sample, unit):
    return has_only_sra_accession(sample, unit) and not config["single_end"]["activate"]

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
    return multiqc_input

# Original from RNA-seq

def get_map_reads_input_R1(wildcards):
    if not is_activated("mergeReads"):
        if config["trimming"]["activate"]:
            if is_paired_end(wildcards.sample):
                return expand(
                    "results/trimmed/{sample}_{unit}_R1.fastq.gz",
                    unit=units.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )
            else:
                return expand(
                    "results/trimmed/{sample}_{unit}_single.fastq.gz",
                    unit=units.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )                
        unit = units.loc[wildcards.sample]
        if has_only_sra_accession(wildcards.sample, wildcards.unit):
            accession = unit["sra"]
            if not config["single_end"]["activate"]:
                return expand("sra-pe-reads/{accession}_1.fastq.gz", accession=accession)
            if config["single_end"]["activate"]:
                return expand("sra-se-reads/{accession}.fastq.gz", accession=accession)
        sample_units = units.loc[wildcards.sample]
        return sample_units["fq1"]
    unit=units.loc[wildcards.sample]
    if all(pd.isna(unit["fq1"])):
        if not config["single_end"]["activate"]:
            return "results/merged/{sample}_R1.fastq.gz"
        if config["single_end"]["activate"]:
            return "results/merged/{sample}_single.fastq.gz"
    ext = units.loc[wildcards.sample]["fq1"][0]
    if ext.endswith("gz") or config['trimming']['activate']:
        ending = ".gz"
    else:
        ending = ""
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R1.fastq" + ending
    return "results/merged/{sample}_single.fastq" + ending


def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        if not is_activated("mergeReads"):
            if config["trimming"]["activate"]:
                return expand(
                    "results/trimmed/{sample}_{unit}_R1.fastq.gz",
                    unit=units.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )
            unit = units.loc[wildcards.sample]
            if has_only_sra_accession(wildcards.sample, wildcards.unit):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                if not config["single_end"]["activate"]:
                    return expand("sra-pe-reads/{accession}_2.fastq.gz", accession=accession)
            sample_units = units.loc[wildcards.sample]
            return sample_units["fq2"]
        unit=units.loc[wildcards.sample]
        if all(pd.isna(unit["fq1"])):
            if not config["single_end"]["activate"]:
                return "results/merged/{sample}_R2.fastq.gz"
            if config["single_end"]["activate"]:
                return ""            
        ext = units.loc[wildcards.sample]["fq1"][0]
        if ext.endswith("gz") or config['trimming']['activate']:
            ending = ".gz"
        else:
            ending = ""
        return ("results/merged/{sample}_R2.fastq" + ending,)
    return ""


def get_star_output_all_units(wildcards, fi,orig=False):
    if fi == "bam":
        resfold = "star"
        outfile = "Aligned.sortedByCoord.out.bam"
    elif fi == "SJ":
        resfold = "star"
        outfile = "SJ.out.tab"
    elif fi == "rsem":
        resfold = "rsem"
        outfile = "mapped.genes.results"
    elif fi == "log":
        resfold = "star"
        outfile = "Log.final.out"
    res = []
    for unit in units.itertuples():
        if is_paired_end(unit.sample_name):
            lib = "pe"
        else:
            lib = "se"
        if is_activated("2nd_pass") and not orig:
            lib=lib+"2"
        if is_activated("mergeReads"):
            res.append(
                "results/{}/{}/{}/{}".format(
                    resfold, lib, unit.sample_name, outfile
                )
            )
        else:
            res.append(
                "results/{}/{}/{}-{}/{}".format(
                    resfold, lib, unit.sample_name, unit.unit_name, outfile
                )
            )
    return list(set(res))


def get_star_bam_uns(wildcards, original=False):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass") and original == False:
        lib = lib + "2"
    if is_activated("mergeReads"):
        return "results/star/{}/{}/Aligned.out.bam".format(
            lib, wildcards.sample
        )
    else:
        return "results/star/{}/{}-{}/Aligned.out.bam".format(
            lib, wildcards.sample, wildcards.unit
        )


def get_star_bam(wildcards, original=False):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass") and not original:
        lib = lib + "2"
    if is_activated("mergeReads"):
        return "results/star/{}/{}/Aligned.sortedByCoord.out.bam".format(
            lib, wildcards.sample
        )
    else:
        return "results/star/{}/{}-{}/Aligned.sortedByCoord.out.bam".format(
            lib, wildcards.sample, wildcards.unit
        )

def get_star_bam_bai(wildcards, original=False):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass") and original == False:
        lib = lib + "2"
    if is_activated("mergeReads"):
        return "results/star/{}/{}/Aligned.sortedByCoord.out.bam.bai".format(
            lib, wildcards.sample
        )
    else:
        return "results/star/{}/{}-{}/Aligned.sortedByCoord.out.bam.bai".format(
            lib, wildcards.sample, wildcards.unit
        )

def get_star_transcript_bam(wildcards):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass"):
        lib = lib + "2"
    if is_activated("mergeReads"):
        return "results/star/{}/{}/Aligned.toTranscriptome.out.bam".format(
            lib, wildcards.sample
        )
    else:
        return "results/star/{}/{}-{}/Aligned.toTranscriptome.out.bam".format(
            lib, wildcards.sample, wildcards.unit
        )

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


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_fastqs_gz(wc):
    if config["trimming"]["activate"]:
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
        if not config["single_end"]["activate"]:
            return expand(
                "sra-pe-reads/{accession}_{read}.fastq.gz", accession=accession, read=wc.read[-1]
            )
        if config["single_end"]["activate"]:
            return expand(
                "sra-se-reads/{accession}.fastq.gz", accession=accession
            )
    fq = "fq{}".format(wc.read[-1])
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
    return config["diffexp"]["contrasts"][wildcards.contrast]

#Returns a path with {sample}-{unit} if merge is deactivated, and only {sample} if its activated in place for {}
def path_merged_cond(path):
    if is_activated("mergeReads"):
        path = path.replace('?','{sample}')
    else:
        path = path.replace('?','{sample}-{unit}')
    if is_activated("2nd_pass"):
        path = path.replace('¿',"2")
    else:
        path = path.replace('¿',"")
    return path

def path_merged_cond_mqc(path):
    if is_activated("mergeReads"):
        return path.replace('?','{unit.sample_name}')
    else:
        return path.replace('?','{unit.sample_name}-{unit.unit_name}')

def get_bg_str1(wc):
    if get_strandedness(units)[0] == 0.5:
        return path_merged_cond("results/bw_uns/?/Signal.Unique.str1.out.bg.sorted")
    else:
        return path_merged_cond("results/bw_str/?/Signal.Unique.str1.out.bg.sorted")

# def get_single_input(wildcards):
#     resfold = "rsem"
#     outfile = "mapped.genes.results"
#     if is_paired_end(wildcards.sample):
#         lib = "pe"
#     else:
#         lib = "se"
#     if is_activated("2nd_pass"):
#         lib=lib+"2"
#     return "results/{}/{}/{}/{}".format(
#             resfold, lib, wildcards.sample, outfile
#         )

def get_single_input(wildcards):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass"):
        lib = lib + "2"
    return "results/rsem/{}/{}/mapped.genes.results".format(
        lib, wildcards.sample
    )

def get_star_log(wildcards):
    if is_paired_end(wildcards.sample):
        lib = "pe"
    else:
        lib = "se"
    if is_activated("2nd_pass"):
        lib = lib + "2"
    return "results/star/{}/{}/Log.final.out".format(
        lib, wildcards.sample
    )

def get_deseq2_end(wc):
    for unit in units.itertuples():
        if is_paired_end(unit.sample_name):
            return "FALSE"
        return "TRUE"