# Simulation pipeline

Simulation pipeline for ecDNA structures.<br/>
The pipeline generates:
- ecDNA template `.fasta` based on a bed file definition
- reads based on the ecDNA template `.fastq`

With the raw data is performed:
- assembly
- evaluation of the assembly
- mapped, sv, cnv calling

## Initialize index and dictionary
If you clone this repo, index and dictionary are already included. If you use different genome or change location of genome please run initialize rule first

```commandline
snakemake --use-conda -s rules/genome_index.smk 
```
## Run entire pipeline using snakemake per sample

```commandline
snakemake --cores 8 \
--use-conda \
--configfile configs/config.yaml \
--config input=data/raw/AnBC.bed name=AnBC outputdir=data/process/AnBC
```

## Simulate reads for ecDNA templates using slurm

For single ecDNA template:

```commandline
#!/bin/bash
#
#SBATCH --job-name=snakemake_main_job
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --time=14-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --output=slurm_logs/%x-%j.log

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date

snakemake --cores 32 \
--use-conda \
--configfile configs/config.yaml \
--config input=data/raw/AnBC.bed name=AnBC outputdir=data/process/AnBC runmode=simulate-mapping-sv

date
```


For multiple templates at once:

```commandline
#!/bin/bash
#
#SBATCH --job-name=snakemake_main_job
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --time=14-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --output=slurm_logs/%x-%j.log

mkdir -p slurm_logs
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

date

snakemake --cores 32 \
--use-conda \
--configfile configs/config.yaml \
--config batch=<dir_batch>  outputdir=data/process/<batch_number> runmode=simulate-mapping-sv

date
```
Argument `batch=<dir_batch>` is set to run the pipeline for multiple templates at once, located under `<dir_batch>`.

## Useful commands

Plot pipeline dag:

```commandline
snakemake --dag | dot -Tpdf > dag.pdf
snakemake --rulegraph | dot -Tpdf > dag_simplified.pdf
```

## Install Quast

```commandline
conda install quast

# GRIDSS (needed for structural variants detection)
quast-download-gridss

# install gene annotation
quast-download-busco

```

## Generate kmers for winnowmap2

Precompute the high frequency k-mers for the different assemblies.

```
meryl count k=15 output merylDB ref.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
```
