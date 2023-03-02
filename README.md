# Snakemake workflow: RNA VariantCalling Docker
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.1-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow performs variant calling and annotation for bulk RNA-Seq samples.

## Authors

* [Matteo Massidda](https://github.com/massiddamt), University of Sassari

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=massiddamt/rna_vc_docker).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## INSTRUCTIONS
Create a virtual environment with the command:
```commandline
mamba create -c bioconda -c conda-forge --name snakemake snakemake=7.18.1 snakedeploy
```
and activate it:
```commandline
conda activate snakemake
```
You can perform the pipeline deploy defining a directory `my_dest_dir` for analysis output and a pipeline tag for a specific version:
```bash
snakedeploy deploy-workflow https://github.com/massiddamt/rna_vc_docker.git 
                    my_desd_dir 
                    --tag v1.0.1
```
To run the pipeline, go inside the deployed pipeline folder and use the command:
```bash
snakemake --use-conda -p --cores all
```