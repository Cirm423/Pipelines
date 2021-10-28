#!/bin/bash

# Set up sbatch parameters

#SBATCH -J snake-pipe
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p q1

# Activating conda and the necessary env

source /data1/share/dcyleung/miniconda3/etc/profile.d/conda.sh

conda activate snakemake

# Making sure Multiqc works (it requires the lang to be utf and not ascii)

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Setting the resources folder for snakemake

export SNAKEMAKE_OUTPUT_CACHE=/data1/share/dcyleung/Pipeline/resources

# Create the directory where the slurm logs from snakemake subtasks will go

mkdir -p logs/slurm

#Running snakemake

#snakemake --cluster 'sbatch -c {threads} -o logs/slurm/{rule}.{wildcards}.out -e logs/slurm/{rule}.{wildcards}.out' -j 10 --latency-wait 60 --use-conda --cache --rerun-incomplete

snakemake --profile /data1/share/dcyleung/Pipeline/snakemake/snake_config
