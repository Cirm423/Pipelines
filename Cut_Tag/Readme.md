# Table of contents

- [Table of contents](#table-of-contents)
- [Starting](#starting)
- [Configuring Snakemake](#configuring-snakemake)
  - [Input dataset](#input-dataset)
  - [Assembly](#assembly)
    - [Spike-in assembly](#spike-in-assembly)
- [Running Snakemake](#running-snakemake)
- [Output](#output)
  - [Plot explanations](#plot-explanations)
    - [Phantompeaks](#phantompeaks)
    - [Deeptools](#deeptools)
- [Too long, don't want to read](#too-long-dont-want-to-read)

# Starting

This readme provides instructions on how to run the snakemake Cut&Tag pipeline on a cluster.

To begin you will need to grab a copy of the pipeline from the shared folder in the cluster:
> hpc2:\
> /home/share/dcyleung/snakemake/Cut_Tag/\
> biocrfhpc1:\
> /home4/share/dcyleung/snakemake/Cut_Tag/\
> biocrfhpc2:\
> /data1/share/dcyleung/Pipeline/snakemake/Cut_Tag/

and place it in your own folder, preferentially in its own folder since snakemake will create many new files.

You will see 2 folders, one named workflow (which includes the code and scripts, please ***do not*** modify anything in this folder) and a config folder; and the run_snk.sh script that will be used to send the pipeline to slurm.

The config folder includes the following files:

- config.yaml
- samples.tsv
- units.tsv
- pe_bamtools_filtering_rules.json
- se_bamtools_filtering_rules.json

Due to the difference in the folder structure in the clusters, the config file of the pipeline is configured for a specific cluster and cannot be used with the other one. Particularly, config.yaml and run_snk.sh are specific to the cluster, while samples.tsv, units.tsv and the files in workflow folder can be used in any cluster.


# Configuring Snakemake


## Input dataset

The previously mentioned *samples.tsv* has to be modified according to your own samples. It is recommended to modify the tsv files in Excel or similar programs since a misplaced tabulation will cause snakemake to improperly reading the config file and probably stop it from working. If you use Excel, please **make sure you save the file as a csv or tsv file** and not as an Excel file, otherwise the pipeline cannot read the file.

The first file you need to modify is *samples.tsv*. It is a tab separated file that looks like this:

> | sample_name | group | batch_effect | control | antibody |
> ------------|---------|--------------|---------|----------|
> | J1_Va_Rep1  | J1 | batch1 |  | IgG  |
> | J1_Va_Rep2	| J1 | batch1 |  | IgG  |
> | J1_Va_Rep3	| J1 | batch1 |  | IgG  |
> | DnmtDKO_Va_Rep1	| DnmtDKO | batch2 | J1_Va_Rep1 | H3K27me3 |
> | DnmtDKO_Va_Rep2	| DnmtDKO | batch2 | J1_Va_Rep2 | H3K27me3 |
> | DnmtDKO_Va_Rep3	| DnmtDKO | batch2 | J1_Va_Rep3 | H3K27me3 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 | J1_Va_Rep1 | H3K27me3 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 | J1_Va_Rep2 | H3K27me3 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 | J1_Va_Rep3 | H3K27me3 |

You need to modify this file to include any samples you want to analyze in the pipeline, along with their group (the condition that will be used in Deseq2 model and for consensus peak calling), batch and control samples. Note that **samples without control will be considered as controls in the pipeline.** Additionally, control samples need to have the antibody as `IgG` (not case sensitive), so that the pipeline can identify and remove these samples from consensus analysis.

**It is also advisable to avoid special characters (like - or |) in the name of the samples as some of them are used by the pipeline to process results, but the pipeline should still work with them.**

**Naming for consensus analysys:** If you want to run consensus peak analysis you need to name your samples in an specific way so the pipe can recognize biological replicates and group them. The way to name them is:

`<sample>_R<replicate>`

where `<sample>` is the the shared name between the biological replicates and `<replicate>` is the specific replicate identifier (usually a number, but can be anything). If you look in the table above, the samples are `J1_Va`, `DnmtDKO_Va` and `DnmtTKO_Va`, each having three replicates `ep1`, `ep2` and `ep3`.

The next file that needs to be modified is *units.tsv*, where you indicate the location of your fastq.gz files. The unit column refer to technical replicates of a sample, e.g. lanes in sequencing. This file looks like this:

> | sample_name |	unit | fragment_len_mean | fragment_len_sd | fq1 | fq2 | sra | platform |
> --------------|--------|-------------------|-----------------|-----|----|------|----------|
> | J1_Va_Rep1  | 1 | | | J1_1.lane1.R1.fastq.gz | J1_1.lane1.R2.fastq.gz | | ILLUMINA |
> | J1_Va_Rep2  | 1 | | | J1_2.lane1.R1.fastq.gz | J1_2.lane1.R2.fastq.gz | | ILLUMINA |
> | J1_Va_Rep3	| 1 | | | J1_3.lane1.R1.fastq.gz | J1_3.lane1.R2.fastq.gz| | ILLUMINA |
> | DnmtDKO_Va_Rep1	| 1 | | |  |  | SRR10194959 | ILLUMINA |
> | DnmtDKO_Va_Rep2	| 1 | | |  |  | SRR10194960 | ILLUMINA |
> | DnmtDKO_Va_Rep3	| 1 | | |  |  | SRR10194961 | ILLUMINA |
> | DnmtTKO_Va_Rep1	| 1 | | |  |  | SRR10194962 | ILLUMINA |
> | DnmtTKO_Va_Rep2	| 1 | | |  |  | SRR10194963 | ILLUMINA |
> | DnmtTKO_Va_Rep3	| 1 | | |  |  | SRR10194964 | ILLUMINA |

You will need to fill this file with either the location of your fastq.gz files or an sra ID for public samples. The path to your files can be the full path to your files, i.e:

> /root/user/snakemake/samples/sample_1/sample_1_R1.fastq.gz

or a relative path from where snakemake is run, which is the directory where the folders *workflow* and *config* lie, i.e:

> samples/sample_1/sample_1_R1.fastq.gz

This latter approach is the preferred one.

**If only the column fq1 is filled, snakemake will try to run the pipeline as single end. If both fq1 and fq2 are filled, snakemake will run the pipeline as paired end. Whether SRA accession samples are considered paired or single end is determined by the *single_end* setting activation in *config.yaml*. This setting should be set accordingly even for fq only runs as some steps make use of it. If both SRA and fastq.gz are present, snakemake will use the fastq.**

The fragment_len_mean and fragment_len_sd refer to the sequencing fragments mean and standard deviation, they can be put in the pipeline if known but are not necessary and the pipeline still doesn't consider them.

Lastly, the *config.yaml* file sets what analyses snakemake will do. This file has been commented to explain what each setting does, so modify the settings to your needs. 

However, some options require further explanation that can be seen bellow.

## Assembly

To indicate which genome assembly you want to use for the reference genome you have to modify the assembly field in *config.yaml*. Currently the pipeline supports UCSC and Genecode assemblies. To use a Genecode assembly you have to provide the Genecode reference to the assembly field. 

Four Genecode assemblies are supported:


> | Organism |	assembly  
> --------------|--------------
> | Human   | GRCh38
> | Human	| GRCh37
> | Mouse	| GRCm39 
> | Mouse	| GRCm38


If you indicate any assembly name that does not appear in the table, the pipeline will try to download an UCSC reference genome. For example, if you want to use the reference mouse genome from UCSC, you can indicate mm39 instead of GCRm39 in the assembly field. 

For organisms other than human or mouse simply indicate the assembly name present in UCSC browser.

### Spike-in assembly

In order to normalize the samples in Cut&Tag it is needed to map the samples against an *E. coli* genome, as the experimental process introduces *E. coli* DNA fragments into the samples at a constant pace. Therefore, these fragments can be used to normalize and then compare different samples.

The recommended genome assembly for the *E. coli* genome is `ASM584v2`, but you can change this assembly to any other of your preference in `config.yaml`, with the `spike_assembly` setting.

# Running Snakemake

To run the pipeline simply send the run_snk.sh script to slurm by doing:

    sbatch run_snk.sh

from the directory where the folders you copied, config and workflow, are located. If run_snk.sh is run from another directory snakemake will not detect the pipeline and it will not work.

Inside *run_snk.sh* you can modify the partition to which snakemake will send the jobs by changing the line `#SBATCH -p q1` and the option `-p q1` in the `snakemake` line if its present (in hpc2 and biocrfhpc1). snakemake will run on a default profile settings in every cluster, if you want to use the another partition, comment (using #) the snakemake line with the profile and uncomment the longer snakemake line (by removing #). Then set your partition by changing or adding `-p <partition>` after the `sbatch` command present in the line (within the single quotes).

Additionally, snakemake is configured to run up to 5 simultaneous jobs by default, but if you're using the pipeline in a cluster with number of jobs restrictions you may want to reduce it to 2 or 1 jobs by changing the `-j` parameter in the `snakemake` line. Snakemake will take at least 2 job spots, one for the pipeline manager and another for the actual jobs, when `-j` is set to 1. In hpc2, `-j` is set to 1 by default since its the cluster with the most restrictions. Therefore if you use the pipeline in the biocrfhpc clusters it will be faster.

# Output

Snakemake will store all the output files in a directory called results. The outputs are organized by steps in the pipeline and samples, so you will se folders like bwa, fastqc, etc; with sample folders within them. For example:

    results/bwa/sample1/bwa_output_files

In results, the qc folder will contain the files `multiqc_report_data` and `multiqc_report.html`, which includes fastqc and rseqc, along with the folder ATACseqQC, where all its plots are found. You can download these files and view in a browser. The file `report.html` will also be found in the root directory of the pipeline, in which you can find the rest of the QC and stats that are not included in multiqc.

The peaks and related info will be in a folder called SEACR (the peak caller program), and all the plots will be included in the pipeline report, which can be found in the root directory of the pipeline when the pipeline is done.

Additionally logs for each step will be stored in the logs folder. 

## Plot explanations

Most figures and scores have explanations within their files on how to interpret them, here are some graphs for additional help.

### Phantompeaks

To get strand cross-correlation scores, as good quality ChIP-seq should produce clusters of enriched DNA sequence tags around sites of interest.

The important scores are nsc and rsc, both explained in output files. The higher, the better. Very good if both are > 1. These scores are present in the 3rd and 2nd from the last columns in *sample*.phantompeak.spp.out respectively, as well as in their individual files. The last column is a quality score based on rsc, which ranges between -2 and 2, lower score meaning worse quality.

### Deeptools

For measuring cumulative enrichment, determines how well a ChIP-seq sample can be distinguished from its control sample (meaning that we can separate ChIP signal from the noise).

- Fingerprint plot

![image](Images/fingerprint.png)
From Deeptools manual

# Too long, don't want to read

The basic steps to run the pipe are:

- Make a copy of the pipeline.
- Put your samples and groups in samples.tsv.
- Put your units and file paths/sra codes in units.tsv
- Change config.yaml to single or paired end, set assembly and any specific step option you need.
- Use `sbatch` to send run_snk.sh to the cluster.
