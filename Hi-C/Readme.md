# Table of contents

- [Table of contents](#table-of-contents)
- [Starting](#starting)
- [Configuring Snakemake](#configuring-snakemake)
  - [Input dataset](#input-dataset)
  - [Mode of use](#mode-of-use)
  - [Restriction Enzymes](#restriction-enzymes)
  - [Assembly](#assembly)
  - [Warning](#warning)
- [Running Snakemake](#running-snakemake)
- [Output](#output)
- [Possible errors](#possible-errors)
- [Too long, don't want to read](#too-long-dont-want-to-read)
- [References](#references)

# Starting

This readme provides instructions on how to run the snakemake Hi-C pipeline on a cluster.

To begin you will need to grab a copy of the pipeline from the shared folder in the cluster:
> hpc2:\
> /home/share/dcyleung/snakemake/Hi-C/\
> biocrfhpc1:\
> /home4/share/dcyleung/snakemake/Hi-C/\
> biocrfhpc2:\
> /data1/share/dcyleung/Pipeline/snakemake/Hi-C/

and place it in your own folder, preferentially in its own folder since snakemake will create many new files.

You will see 2 folders, one named workflow (which includes the code and scripts, please ***do not*** modify anything in this folder) and a config folder; and the run_snk.sh script that will be used to send the pipeline to slurm.

The config folder includes the following files:

- config.yaml
- samples.tsv
- units.tsv

Due to the difference in the folder structure in the clusters, the config file of the pipeline is configured for a specific cluster and cannot be used with the other one. Particularly, config.yaml and run_snk.sh are specific to the cluster, while samples.tsv, units.tsv and the files in workflow folder can be used in any cluster.


# Configuring Snakemake


## Input dataset

The previously mentioned *samples.tsv* has to be modified according to your own samples. It is recommended to modify the tsv files in Excel or similar programs since a misplaced tabulation will cause snakemake to improperly reading the config file and probably stop it from working. If you use Excel, please **make sure you save the file as a csv or tsv file** and not as an Excel file, otherwise the pipeline cannot read the file.

The first file you need to modify is *samples.tsv*. It is a tab separated file that looks like this:

> | sample_name | group | batch_effect |
> ------------|---------|--------------|
> | J1_Va_Rep1  | J1 | batch1 |
> | J1_Va_Rep2	| J1 | batch1 |
> | J1_Va_Rep3	| J1 | batch1 |
> | DnmtTKO_Co_Rep1	| DnmtTKO | batch2 |
> | DnmtTKO_Co_Rep2	| DnmtTKO | batch2 |
> | DnmtTKO_Co_Rep3	| DnmtTKO | batch2 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 |
> | DnmtTKO_Va_Rep1	| DnmtTKO | batch3 |

You need to modify this file to include any samples you want to analyze in the pipeline, along with their group, which represents biological replicates of the same sample. **Samples belonging to the same group can be merged to create a Hi-C matrix for all the biological replicates together.** The batch effect column currently does nothing, but it is good to keep track of it.

**It is also advisable to avoid special characters (like - or |) in the name of the samples as some of them are used by the pipeline to process results, but the pipeline should still work with them.**

The next file that needs to be modified is *units.tsv*, where you indicate the location of your fastq.gz files. The unit column refer to technical replicates of a sample, e.g. lanes in sequencing. This file looks like this:

> | sample_name |	unit | fq1 | fq2 | sra | platform |
> --------------|--------|-----|----|------|----------|
> | J1_Va_Rep1  | 1 | J1_1.lane1.R1.fastq.gz | J1_1.lane1.R2.fastq.gz | | ILLUMINA |
> | J1_Va_Rep2  | 1 | J1_2.lane1.R1.fastq.gz | J1_2.lane1.R2.fastq.gz | | ILLUMINA |
> | J1_Va_Rep3	| 1 | J1_3.lane1.R1.fastq.gz | J1_3.lane1.R2.fastq.gz| | ILLUMINA |
> | DnmtTKO_Co_Rep1	| 1 |  |  | SRR10194959 | ILLUMINA |
> | DnmtTKO_Co_Rep2	| 1 |  |  | SRR10194960 | ILLUMINA |
> | DnmtTKO_Co_Rep3	| 1 |  |  | SRR10194961 | ILLUMINA |
> | DnmtTKO_Va_Rep1	| 1 |  |  | SRR10194962 | ILLUMINA |
> | DnmtTKO_Va_Rep2	| 1 |  |  | SRR10194963 | ILLUMINA |
> | DnmtTKO_Va_Rep3	| 1 |  |  | SRR10194964 | ILLUMINA |

You will need to fill this file with either the location of your fastq.gz files or an sra ID for public samples. The path to your files can be the full path to your files, i.e:

> /root/user/snakemake/samples/sample_1/sample_1_R1.fastq.gz

or a relative path from where snakemake is run, which is the directory where the folders *workflow* and *config* lie, i.e:

> samples/sample_1/sample_1_R1.fastq.gz

This latter approach is the preferred one.

**Note that this pipeline currently only supports Paired End samples as most methods currently available only accept Paired End samples. Running it will Single End samples will probably crash the pipeline**

Lastly, the *config.yaml* file sets what analyses snakemake will do. This file has been commented to explain what each setting does, so modify the settings to your needs.

Hi-C many different parameters for each step, each one of them is briefly explained in the config.yaml file itself. As this pipeline uses [fanC](https://github.com/vaquerizaslab/fanc)[1] as a base, it is recommended that you also check their [documentation](https://fan-c.readthedocs.io/en/latest/)

## Mode of use

This pipeline can be used in different ways. The main "switches" are found in the fanc section of config.yaml: `merge_groups` at the start of the fanc section and `activate` in the analysis section.

- If you only want to get the `.hic` files and do your own analysis and plots, turn off the analysis activation so the pipeline only maps your fastqs and generates the matrix, either for merged groups or single samples.

- If you want to analyze your data is advised that you first check how similar your biological replicates are in order to be merged safely. To do this, run the pipeline with merge_groups as `False` with the analysis activated. In this case, each sample will have its own Hi-C matrix generated and compared by doing a PCA plot.

- If you are confident in your biological replicates or the PCA looks ok, you can then run the analysis with merge_groups as `True` and the analysis activated, which will perform the rest of the analysis.

## Restriction Enzymes

In Hi-C, the genome needs to be "digested" and separated into fragments that change depending on the restriction enzymes utilized by the protocol
used to sequence your samples. The pipeline can digest the genome if you provide the name of the restriction enzymes as long as they're found in [REBASE](http://rebase.neb.com/rebase/rebase.html). By default the pipeline is set to digest using the restriction enzymes used by the common Arima protocol.

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

## Warning

This pipeline can take a long time to run (a week or more) depending on your settings.

**Using the full genome can make the Hi-C matrix generation step really slow, so if you are interested only in specific chromosomes or regions, it is better to restrict the pipeline to only construct Hi-C matrices of the desired regions in the `chr` option of the fanc section.**

**Setting a high resolution matrix (small bp bins) for the analysis can make some steps of the analysis very slow. If you want to test your data and need some quick plots, maybe is better to choose a lower resolution.**

**This pipeline has many parameters for every step, and Hi-C can be very sensitive to them, so it is recommended that you read the documentation for every step or the explanations that are in `config.yaml` at the very least.**

# Running Snakemake

To run the pipeline simply send the run_snk.sh script to slurm by doing:

    sbatch run_snk.sh

from the directory where the folders you copied, config and workflow, are located. If run_snk.sh is run from another directory snakemake will not detect the pipeline and it will not work.

Inside *run_snk.sh* you can modify the partition to which snakemake will send the jobs by changing the line `#SBATCH -p q1` and the option `-p q1` in the `snakemake` line if its present (in hpc2 and biocrfhpc1). snakemake will run on a default profile settings in every cluster, if you want to use the another partition, comment (using #) the snakemake line with the profile and uncomment the longer snakemake line (by removing #). Then set your partition by changing or adding `-p <partition>` after the `sbatch` command present in the line (within the single quotes).

Additionally, snakemake is configured to run up to 5 simultaneous jobs by default, but if you're using the pipeline in a cluster with number of jobs restrictions you may want to reduce it to 2 or 1 jobs by changing the `-j` parameter in the `snakemake` line. Snakemake will take at least 2 job spots, one for the pipeline manager and another for the actual jobs, when `-j` is set to 1. In hpc2, `-j` is set to 1 by default since its the cluster with the most restrictions. Therefore if you use the pipeline in the biocrfhpc clusters it will be faster.

# Output

Snakemake will store all the output files in a directory called results. The outputs are organized by steps in the pipeline and samples, so you will se folders like bwa, fastqc, etc; with sample folders within them. For example:

    results/bwa/sample1_bwa_outputs

In results, the qc folder will contain the files `multiqc_report_data` and `multiqc_report.html`, which includes fastqc and rseqc, along with the folder ATACseqQC, where all its plots are found. You can download these files and view in a browser. The file `report.html` will also be found in the root directory of the pipeline, in which you can find the rest of the QC and stats that are not included in multiqc.

In addition to the plots and qc that will be present in these files, the pipeline will provide you with 3 Hi-C matrices files:

  - A FanC Hi-C file that can be used with FanC software to make plots.
  - A Juicer Hi-C file with multiple resolutions that can be used with the juicebox software tools.
  - A cooler Hi-C file that can be used with cooler and domaincaller.

Additionally logs for each step will be stored in the logs folder. 

# Possible errors

The steps fanc_compartments, fanc_expected and fanc_loops_annotate are a bit sensible and can fail since all of them could be accessing the Hi-C matrix at the same time, especially if they use the same node as they can have memory issues. The pipeline tries to go over this by repeating the step in case it fails. If despite this the pipeline fails to complete with an error in these steps just resubmitting the pipeline to slurm fixes the problem in most cases. **Check that the pipeline has finished running all current steps before resubmitting.**


# Too long, don't want to read

The basic steps to run the pipe are:

- Make a copy of the pipeline.
- Put your samples and groups in samples.tsv.
- Put your units and file paths/sra codes in units.tsv
- Change config.yaml to your restriction enzymes, set assembly and any specific step option you need.
- Use `sbatch` to send run_snk.sh to the cluster.

# References
[1] Kruse, K., Hug, C.B. & Vaquerizas, J.M. 
FAN-C: a feature-rich framework for the analysis and visualisation of chromosome conformation capture data. 
Genome Biol 21, 303 (2020). 
https://doi.org/10.1186/s13059-020-02215-9