# Table of contents

- [Table of contents](#table-of-contents)
- [Starting](#starting)
- [Configuring Snakemake](#configuring-snakemake)
  - [Input dataset](#input-dataset)
  - [Assembly](#assembly)
  - [Deseq2 Model](#deseq2-model)
  - [Transposable Elements (TE)](#transposable-elements-te)
  - [Single Sample](#single-sample)
    - [Genes](#genes)
    - [TEs](#tes)
- [Running Snakemake](#running-snakemake)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
  - [Pipeline not starting](#pipeline-not-starting)
  - [Some step in the pipeline fails](#some-step-in-the-pipeline-fails)
- [Acknowledgements](#acknowledgements)

# Starting

This readme provides instructions on how to run the snakemake RNA-seq pipeline on a cluster.

To begin you will need to grab a copy of the pipeline from the shared folder in the cluster:
> hpc2:\
> /home/share/dcyleung/snakemake/RNA-seq/\
> biocrfhpc1:\
> /home4/share/dcyleung/snakemake/RNA-seq/\
> biocrfhpc2:\
> /data1/share/dcyleung/Pipeline/snakemake/RNA-seq/

and place it in your own folder, preferentially in its own folder since snakemake will create many new files.

You will see 2 folders, one named workflow (which includes the code and scripts, please ***do not*** modify anything in this folder) and a config folder; and the run_snk.sh script that will be used to send the pipeline to slurm.

The config folder includes the following files:

- config.yaml
- samples.tsv
- units.tsv

Due to the difference in the folder structure in the clusters, the config file of the pipeline is configured for a specific cluster and cannot be used with the other one. Specifically, config.yaml and run_snk.sh are specific to the cluster, while samples.tsv, units.tsv and the files in workflow folder can be used in any cluster.

Before starting, make sure that your fastq files are compressed in fastq.gz format, both to save space in clusters and to make sure the pipeline works properly.

# Configuring Snakemake


## Input dataset

The previously mentioned *samples.tsv* has to be modified according to your own samples. It is recommended to modify the tsv files in Excel or similar programs since a misplaced tabulation will cause snakemake to improperly reading the config file and probably stop it from working. If you use Excel, please **make sure you save the file as a csv or tsv file** and not as an Excel file, otherwise the pipeline cannot read the file.

The first file you need to modify is *samples.tsv*. It is a tab separated file that looks like this:

> | sample_name |	condition |
> ------------|-----------------
> | A1  | treated |
> | B1	| untreated |
> | A2	| treated |
> | B2	| untreated |

You need to modify this file to include any samples you want to analyze in the pipeline, and the condition that will be used in Deseq2 model. If you are not going to use Deseq2 leave condition as *treated* for every sample, but the samples should still be filled in this file.

If you need to analyze more conditions in Deseq2, you can add more columns to the right of *condition* column and set each sample treatment in that column.


The next file that needs to be modified is *units.tsv*, where you indicate the location of your fastq.gz files. The unit_name columns refer to technical replicates of a sample, e.g. lanes in sequencing. Different units from the same sample will be merged together (eg. Sample A1 lane1 and lane2 fq files will be merged). If your samples don't have technical replicates, give the same unit_name to all samples. This file looks like this:

> | sample_name |	unit_name | fq1 | fq2 | sra | adapters | strandedness |
> ------------|---------------|-----|-----|-----|----------|--------------|
> | A1  | lane1 | A1.lane1.R1.fastq.gz | A1.lane1.R2.fastq.gz | | | 1 |
> | A1  | lane2 | A1.lane2.R1.fastq.gz | A1.lane2.R2.fastq.gz | | | 1 |
> | B1	| lane1 | B1.lane1.R1... | B1.lane1.R2...| | | 1 |
> | B1	| lane2 | ... | ... | | | 1 |
> | A2	| lane1 | ... | ... | | | 1 |
> | A2	| lane2 | ... | ... | | | 1 |
> | B2	| lane1 | ... | ... | | | 1 |
> | B2	| lane2 | ... | ... | | | 1 |

You will need to fill this file with either the location of your fastq.gz files or an sra ID for public samples. The path to your files can be the full path to your files, i.e:

> /root/user/snakemake/samples/sample_1/sample_1_R1.fastq.gz

or a relative path from where snakemake is run, which is the directory where the folders *workflow* and *config* lie, i.e:

> samples/sample_1/sample_1_R1.fastq.gz

This last approach is the preferred one.

**If only the column fq1 is filled, snakemake will run the pipeline as single end. If both fq1 and fq2 are filled, snakemake will run the pipeline as paired end. Whether SRA accession samples are considered paired or single end is determined by the *single_end* setting activation in *config.yaml*. This setting has no effect on local samples. If both SRA and fastq.gz are present, snakemake will use the fastq.**

The last necessary column is strandedness, which **needs to be equal for all units**. The number equivalences can be found here:  

> Number | Strandedness |
> -------|---------------
> 1 | Stranded |
> 0.5 | Unstranded |
> 0 | Reverse Stranded* |
> *Only used in Illumina TruSeq protocol  

             
The adapters column is only used for trimming and can be left empty. However if you need to trim an adapter you have to input the sequence of  the adapter in this column for each unit.

Lastly, the *config.yaml* file sets what analyses snakemake will do. This file has been commented to briefly explain what each setting does, so modify the settings to your needs. 

However, some options require further explanation that can be seen bellow.

## Assembly

To indicate which genome assembly you want to use for the reference genome you have to modify the assembly field in *config.yaml*. Currently the pipeline supports UCSC and Genecode assemblies. To use a Genecode assembly you have to provide the Genecode reference to the assembly field. 

Four Genecode assemblies are supported:


> | Organism |	assembly   | Ucsc rmsk
> --------------|---------------|---------
> | Human   | GRCh38 | hg38
> | Human	| GRCh37 | hg19
> | Mouse	| GRCm39 | mm39
> | Mouse	| GRCm38 | mm10

In case you are using the pipeline for TEs, the remasker file indicated in the 3rd column will be downloaded from UCSC to be used as a reference.

If you indicate any assembly name that does not appear in the table, the pipeline will try to download an UCSC reference genome. For example, if you want to use the reference mouse genome from UCSC, you can indicate mm39 instead of GCRm39 in the assembly field. 

For organisms other than human or mouse simply indicate the assembly name present in UCSC browser.

## Deseq2 Model

To use Deseq2 there needs to be a model set in *config.yaml*. The model needs to be set even if you are not going to activate Deseq2, in this case just leave the default model as is. The default model used is 

    ~condition

which uses the condition column present in *samples.tsv*. If you added additional columns to be used for testing in Deseq2, they need to be put in the model. For example, if you added the column *mutated* to *samples.tsv*, the model present in *config.yaml* would be:

    ~condition + mutated

Then Deseq2 will produce a result file by comparing 2 contrasts, which are indicated in *config.yaml* contrasts section. A base one is already provided in *config.yaml*, but more can be added in the following way:

    contrasts:
        treated-vs-untreated:
          - treated
          - untreated
        mutated-vs-control:
          - mutated
          - control    

where mutated-vs-control is what will be used for the file name, and the conditions inside is the contrasts that will be compared in Deseq2. Note that all the contrasts here have to appear in *samples.tsv* for Deseq2 to work.

## Transposable Elements (TE)

Within the DEseq2 config block is also the TE configuration. If you activate the TE option snakemake will produce DEseq2 output for TEs in addition to the normal DEseq2 output using the RepeatMasker file from UCSC for your selected assembly, including the R object and the pca (if activated). Both DEseq2 and TEs will use the same model specified in this config block.

The other setting in this section is the filter, which will filter out elements in DEseq2 with less than the filter number of counts. Is set to 10 by default but you can change it to suit your needs.

## Single Sample

### Genes

To be used when you have no replicates for a sample. It is activated by turning the single option in config.yaml to True, and it will perform an alternative analysis to Deseq2. Therefore, the single option cannot be used at the same time at Deseq2 (and therefore pca), otherwise the pipeline will not work. Single mode also requires mergeReads to be activated, as it will be comparing samples and not units.

Additionally, the treatment of one of the samples in samples.tsv **must be called `control`**. The other treatments do not matter, as every sample will be compared to the control sample.

Once you run the pipeline in single mode with one sample as control, you can rerun it with another sample as control and new files comparing samples to the new control will be generated without needing to run the whole pipeline.

### TEs

In a similar way, you can obtain TE differential expression for samples with no replicates by turning *TE_single* to activate in the `config.yaml` file. TE_single shares the same requirements as single mode does, and both single and TE_single can be activated at the same time, producing output for both genes and TEs.

The only difference for TE_single is the filter, which will remove from the analysis elements with less than the filter number of counts, similarly to the TE analysis with DEseq2.

# Running Snakemake

To run the pipeline simply send the run_snk.sh script to slurm by doing:

    sbatch run_snk.sh

from the directory where the folders you copied, config and workflow, are located. If run_snk.sh is run from another directory snakemake will not detect the pipeline and it will not work. This will run snakemake with the default config for sbatch for the cluster, using the partition general for hpc2 and q1 for biocrf clusters.

Inside *run_snk.sh* you can modify the partition to which snakemake will run by changing the line `#SBATCH -p q1`, which is only one job. As for the partition where snakemake will send the jobs, you can change the last 2 lines of the `run_snk.sh` file. If you want to use the another partition, comment (using #) the snakemake line with the profile and uncomment the longer snakemake line (by removing #). Then set your partition by changing or adding `-p <partition>` after the `sbatch` command present in the line (within the single quotes).

Additionally, snakemake is configured to run up to 10 simultaneous jobs by default in biocrfhpcs, 2 for hpc2; but if you're using the pipeline in a cluster with number of jobs restrictions you may want to reduce it to 2 or 1 jobs by changing the `-j` parameter in the long `snakemake` line. Snakemake will take at least 2 job spots, one for the pipeline manager and another for the actual jobs, when `-j` is set to 1. In hpc2, `-j` is set to 2 by default since its the cluster with the most restrictions. Therefore if you use the pipeline in the biocrfhpc clusters it will be faster.

# Output

Snakemake will store all the output files in a directory called results. The outputs are organized by steps in the pipeline and samples, so you will se folders like star, rsem, etc; with sample folders within them. For example:

    results/rsem/pe/sample1/star_output_files

In results, the qc folder will contain the files `multiqc_report_data` and `multiqc_report.html`, which includes fastqc and rseqc. You can download these files and view in a browser. 

If you activated deseq2, in the deseq2 results folder you can find `all.rds` which contains a Deseq object called dds of all your data. You can load it in R in case you want to further explore the data with Deseq2. The pipeline will also produce some plots like the pca (if activated) so you can initially asses your data.

For single and TE_single modes, the pipeline will create the single and TE_single folders in results respectively, where files named {sample}\_vs_{control}.tsv will be placed, comparing all the samples in `samples.tsv` against the sample marked as control. These files contain the differential expression between each sample and the control. You can rerun the pipeline with a different sample marked as control, which will run only the last step to compare samples against a new control and therefore will be much quicker.

Additionally logs for each step will be stored in the logs folder. 

# Troubleshooting

## Pipeline not starting

Most of these issues have to do with the formatting of units.tsv or samples.tsv, but you may have also activated incompatible config settings.

- 'utf-8' codec can't decode...

If this line appears in the error message in the pipeline slurm output file, it means that there is a special character somewhere in either `units.tsv` or `samples.tsv` (usually the message specifies which file). Libreoffice's calc and nano can show special characters, so you can use them to check your files and remove them. 

Another possibility is that you saved any of the tsv files as an Excel .xml file, in which case they won't work unless you save them as csv files with tab separators.

- Wrong samples or units formatting

If snakemake complains that it cannot find a column in one of the tsv files, or that the values in the column are the wrong type of data; it usually means that there are wrongly placed tabs in the file. You may want to check that the columns are properly formatted in a spreadsheet program and also try to check if there are extra tabs anywhere in the file (like another line at the bottom, or after the last column).

- Directory locked
  
If snakemake complains that the directory is locked, it means that either snakemake is already running or that the previous snakemake job crashed or was cancelled and could not unlock the directory. If snakemake is still running wait until it finishes. If the directory is locked and snakemake is not running, run the following code to unlock the directory (from the RNA-seq directory):

     source ~/../share/dcyleung/miniconda3/etc/profile.d/conda.sh 
     conda activate snakemake 
     snakemake --unlock 
     conda deactivate 
     source ~/.bashrc

## Some step in the pipeline fails

The only step that I've found that can fail that is not my fault (ehem) is rsem. There seems to be a bug where rsem is failing to calculate the confidence intervals (ci), seems to happen more with sra samples. If this is the case, deactivate ci for rsem in `config.yaml`, there are further instructions in the file.

Another step that can fail with a not useful log is Deseq2 init (thanks R), which usually has to do with problems with the model. You can try to run with only `~condition` as a model to check that it works.

If any other step is failing for you, please tell me so I can fix it.

# Acknowledgements

This pipeline is based on the pipeline made by jafors:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5245549.svg)](https://doi.org/10.5281/zenodo.5245549) [Snakemake workflow: rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)