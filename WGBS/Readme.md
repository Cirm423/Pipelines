# Table of contents

- [Table of contents](#table-of-contents)
- [Starting](#starting)
- [Configuring Snakemake](#configuring-snakemake)
  - [Input dataset](#input-dataset)
  - [Assembly](#assembly)
  - [Mode](#mode)
    - [Bismark](#bismark)
    - [Bwa-meth](#bwa-meth)
    - [Which mode to choose](#which-mode-to-choose)
  - [Differential Methylation Analysis](#differential-methylation-analysis)
    - [Co-variates analysis](#co-variates-analysis)
- [Running Snakemake](#running-snakemake)
- [Output](#output)

# Starting

This readme provides instructions on how to run the snakemake WGBS pipeline on a cluster.

To begin you will need to grab a copy of the pipeline from the shared folder in the cluster:
> hpc2:\
> /home/share/dcyleung/snakemake/WGBS/\
> biocrfhpc1:\
> /home4/share/dcyleung/snakemake/WGBS/\
> biocrfhpc2:\
> /data1/share/dcyleung/Pipeline/snakemake/WGBS/

and place it in your own folder, preferentially in its own folder since snakemake will create many new files.

You will see 2 folders, one named workflow (which includes the code and scripts, please ***do not*** modify anything in this folder) and a config folder; and the run_snk.sh script that will be used to send the pipeline to slurm.

The config folder includes the following files:

- config.yaml
- samples.tsv
- units.tsv

Due to the difference in the folder structure in the clusters, the config file of the pipeline is configured for a specific cluster and cannot be used with the other one. Specifically, config.yaml and run_snk.sh are specific to the cluster, while samples.tsv, units.tsv and the files in workflow folder can be used in any cluster.

# Configuring Snakemake


## Input dataset

The previously mentioned *samples.tsv* has to be modified according to your own samples. It is recommended to modify the tsv files in Excel or similar programs since a misplaced tabulation will cause snakemake to improperly reading the config file and probably stop it from working. If you use Excel, please **make sure you save the file as a csv or tsv file** and not as an Excel file, otherwise the pipeline cannot read the file.

The first file you need to modify is *samples.tsv*. It is a tab separated file that looks like this:

> | sample | group | 
> ------------|---------|
> | A1  | treatment |
> | B1	| untreated |
> | A2	| treatment |
> | B2	| untreated |

You need to modify this file to include any samples you want to analyze in the pipeline, along with their group (the condition that will be used in Methylkit).

**It is also advisable to avoid special characters (like - or _) in the name of the samples as some of them are used by the pipeline to process results, but the pipeline should still work with them.**

The next file that needs to be modified is *units.tsv*, where you indicate the location of your fastq.gz files. The unit column refer to technical replicates of a sample, e.g. lanes in sequencing. This file looks like this:

> | sample_name |	unit | fq1 | fq2 | sra |
> --------------|--------|-----|----|------|
> | A1  | lane1 | A1.lane1.R1.fastq.gz | A1.lane1.R2.fastq.gz |  |
> | A1  | lane2 | A1.lane2.R1.fastq.gz | A1.lane2.R2.fastq.gz |  |
> | B1	| lane1 | B1.lane1.R1... | B1.lane1.R2...|  |
> | B1	| lane2 | ... | ... | ... |
> | A2	| lane1 | ... | ... | ... |
> | A2	| lane2 | ... | ... | ... |
> | B2	| lane1 | ... | ... | ... |
> | B2	| lane2 | ... | ... | ... |

You will need to fill this file with either the location of your fastq.gz files or an sra ID for public samples. The path to your files can be the full path to your files, i.e:

> /root/user/snakemake/samples/sample_1/sample_1_R1.fastq.gz

or a relative path from where snakemake is run, which is the directory where the folders *workflow* and *config* lie, i.e:

> samples/sample_1/sample_1_R1.fastq.gz

This last approach is the preferred one.

**If only the column fq1 is filled, snakemake will run the pipeline as single end. If both fq1 and fq2 are filled, snakemake will run the pipeline as paired end. Whether SRA accession samples are considered paired or single end is determined by the *single_end* setting activation in *config.yaml*. If both SRA and fastq.gz are present, snakemake will use the fastq.**

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

## Mode

This pipeline can be run with 2 alternative mappers and methylation callers, and the mode parameter lets you choose which version you want to use. The options are:

### Bismark

[Bismark](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html) is the golden standard for methylation analysis, and will manage everything from mapping to calling methylation to QC. It provides high accuracy with a decent speed.

### Bwa-meth

[Bwa-meth](https://github.com/brentp/bwa-meth) is a newer methylation mapper that provides similar accuracy to bismark but it is much quicker. It uses [MethylDackel](https://github.com/dpryan79/MethylDackel) as it's methylation caller and methylation QC tool. 

### Which mode to choose

If you need to get results as soon as possible, bwa-meth will finish several times faster than bismark. Usually, a bwa-meth run takes 1 to 2 days, while a bismark run can take 4+ days (depending on your sample depth). Both provide similar accuracy so your results should be similar between the two. 
If you want a more standardized method for publication, then you can choose to use the bismark mode instead.
Both modes allow you to use methylkit to detect Diferentially Methylated Regions (DRMs).

## Differential Methylation Analysis

This pipeline can use [MethylKit](https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) to find DRMs. As stated before, output from both Bismark and Bwa-meth can be used with MethylKit.

**If you intend to use bwa-meth with Methylkit, the parameter `methyl_kit` in methyldackel options in `config.yaml` must be set to True.**

Otherwise Methylkit will not run as methyldackel output format will not be able to be read by MethylKit. Bismark does not need any special parameter to be used by MethylKit.

At the moment only CpGs DRMs will be analyzed with MethylKit.

### Co-variates analysis

Methylkit allows to account for co-variate effects to DRMs. If you want to include them in your analysis, you can simply add more columns to `samples.tsv`, which will set Methylkit to run with co-variates and use the new columns as co-variates. For example:

> | sample | group | smoker | Diabetic | Age |
> ------------|---------|------|----|-----|
> | A1  | treatment | No  | No | 32 |
> | B1	| untreated | No  | No | 49 |
> | A2	| treatment | Yes | Yes | 25 |
> | B2	| untreated | Yes | No | 37 |

This table will add the co-variates named "smoker", "Diabetic", and "Age" to the co-variates list in Methylkit, with the values present in their respective column, and run Methylkit taking in account these variables.

Note: This part has not been thoroughly tested so expect some errors, possibly when using numeric variables.

# Running Snakemake

To run the pipeline simply send the run_snk.sh script to slurm by doing:

    sbatch run_snk.sh

from the directory where the folders you copied, config and workflow, are located. If run_snk.sh is run from another directory snakemake will not detect the pipeline and it will not work.

Inside *run_snk.sh* you can modify the partition to which snakemake will send the jobs by changing the line `#SBATCH -p q1` and the option `-p q1` in the `snakemake` line if its present (in hpc2 and biocrfhpc1). snakemake will run on a default profile settings in every cluster, if you want to use the another partition, comment (using #) the snakemake line with the profile and uncomment the longer snakemake line (by removing #). Then set your partition by changing or adding `-p <partition>` after the `sbatch` command present in the line (within the single quotes).

Additionally, snakemake is configured to run up to 5 simultaneous jobs by default, but if you're using the pipeline in a cluster with number of jobs restrictions you may want to reduce it to 2 or 1 jobs by changing the `-j` parameter in the `snakemake` line. Snakemake will take at least 2 job spots, one for the pipeline manager and another for the actual jobs, when `-j` is set to 1. In hpc2, `-j` is set to 1 by default since its the cluster with the most restrictions. Therefore if you use the pipeline in the biocrfhpc clusters it will be faster.

# Output

Snakemake will store all the output files in a directory called results. The outputs are organized by steps in the pipeline, so you will se folders like bismark, methylkit, etc; with sample folders within them. For example:

    results/bismark_mapped/sample1.mapped.bam

In results, the qc folder will contain the files `multiqc_report_data` and `multiqc_report.html`, which includes most of the QC in the pipeline. You can download these files and view in a browser. The report.html created in the root directory of the pipeline can also contain some QC plots depending on the mode set for the run.

You should pay special attention to the Mbias plots, as they can indicate methylation bias that can influence your results. The issues pointed out by the Mbias plots can usually by solved by trimming or ignoring the first and last few BP of the fragments.

In the diffmeth results folder you can find `methDB.RDS` object which contains a methylkit object called methDB with all your files and experimental design already included. You can load it in R in case you want to further explore the data with Methylkit. At the time the methDB object has been stored in the pipeline, no filtering or tiles have been created so you try your own settings. The pipeline will also produce some plots like the pca so you can initially asses your data.

Additionally logs for each step will be stored in the logs folder. 