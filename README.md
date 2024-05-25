# Snakemake RNA VC docker Workflow Execution

This is how to reproduce the workflow from [massiddamt/rna_vc_docker](https://github.com/massiddamt/rna_vc_docker) using SRR390728 and SRR390729.

NCBI metadata:
- [SRR390728](https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR390728&display=metadata)
- [SRR390729](https://www.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR390729&display=metadata)

## Step 1: Install Conda

```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer script and install Miniconda in the home directory
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

# Activate Conda
source $HOME/miniconda/etc/profile.d/conda.sh
```

## Step 2: Install Mamba

```bash
# Use Conda to install Mamba
conda install -y -c conda-forge mamba
```

## Step 3: Install Snakemake and Snakedeploy

```bash
# Create a Conda environment named "snakemake" and install Snakemake and Snakedeploy
mamba create -y -c conda-forge -c bioconda --name snakemake snakemake snakedeploy

# Activate the "snakemake" environment
conda activate snakemake
```

## Step 4: Download data files

```bash
# Install SRA tools to download the files
conda install -c bioconda sra-tools

# Download .sra files
prefetch SRR390728
prefetch SRR390729

# Convert the downloaded .sra files into fastq files
fastq-dump --split-files SRR390728
fastq-dump --split-files SRR390729
```

## Step 5: Prepare input files

1. Change units.tsv file into:

```
sample	unit	fq1	fq2
SRR390728	unit1	/workspace/rna_vc_docker/SRR390728_1.fastq	/workspace/rna_vc_docker/SRR390728_2.fastq
SRR390729	unit2	/workspace/rna_vc_docker/SRR390729_1.fastq	/workspace/rna_vc_docker/SRR390729_2.fastq
```

2. Change samples.tsv file into:

```
sample	odp	units	condition	patient
SRR390728	100	unit1	T	pat_1
SRR390729	100	unit2	T	pat_2
```

## Step 6: Configure workflow

Update paths in the `config.yaml` file:

```yaml
paths:
    workdir: "/workspace/rna_vc_docker"
    results_dir: "/workspace/rna_vc_docker/results"
    tmp_dir: "/workspace/rna_vc_docker/tmp"
```

## Step 7: Run the workflow

```bash
snakemake --cores all --use-conda
```

Finally, the result will be saved at `/workspace/rna_vc_docker/results`.
