# BGCflow
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics)

Snakemake workflow to combine internal &amp; public dataset for downstream Biosynthetic Gene Clusters (BGCs) analyses

## Workflow overview
![dag](workflow/report/images/rulegraph.svg)
## Usage
### Step 1: Clone the workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis. _Make sure to have the right access / SSH Key._

    git clone git@github.com:NBChub/bgcflow.git
    cd bgcflow

### Step 2: Configure workflow
#### Setting Up Your Samples Information
Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.csv` to specify the strains to analyse. Further, update reference `prokka-db.csv` with high quality manually annotated genomes that you would like to use for prokka annoations as priority. 

`samples.csv` example:

| genome_id       | source | organism                        | genus        | species | strain     |
|----------------:|-------:|--------------------------------:|-------------:|--------:| ----------:|
| NBC_01270.1     | custom | Streptomyces sp. NBC_01270      | Streptomyces | sp.     | NBC001270  |
| GCF_000359525.1 | ncbi   | Streptomyces albus strain J1074 | Streptomyces | albus   | J1074      |

`prokka-db.csv` example for Actinomycete group:

| Accession       | Strain Description             |
|----------------:|-------------------------------:|
| GCA_000203835.1 | Streptomyces coelicolor A3(2)  |
| GCA_000196835.1 | Amycolatopsis mediterranei U32 |

Further formatting rules will be defined in the `workflow/schemas/` folder.

### Step 3: Install Snakemake

Installing Snakemake using [Mamba](https://github.com/mamba-org/mamba) is advised. In case you don’t use [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) you can always install [Mamba](https://github.com/mamba-org/mamba) into any other Conda-based Python distribution with:

    conda install -n base -c conda-forge mamba

Then install Snakemake with:

    mamba create -c conda-forge -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via:

    snakemake --use-conda -n

Execute the workflow locally via:

    snakemake --use-conda --cores $N --keep-going

Check you job DAG by executing:

    snakemake --dag | dot -Tsvg > workflow/report/images/dag.svg

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

## General steps for the workflow (To Do)
1. Download public genome fasta files from NCBI using input table public.csv 
2. Download internal genome fasta files from Azure Data Lake using input table internal.csv
3. Run prokka annotation and save output files
4. Run GTDB-tk on internal genomes and download GTDB taxonomy for public genomes
5. Run antiSMASH on prokka generated gbk files
6. (In batch) Run autoMLST-wrapper to get phylogenetic tree
7. (In batch) Run antiSMASH-db import to setup sql server
8. (In batch) Run BiG-SCAPE to get similarity network
9. (In batch) Run BiG-SLICE to get GCFs
10. (In batch) Run BiG-FAM db to setup sql server
11. (In batch) Run Roary to get pangenome presence absence matrix
12. (In batch) Run BGC analytics to get dataframes df_bgcs, df_genomes, df_genes populated with relevant information
13. (In batch) Run BGC analytics to get visualization of pangenome, bgc data
14. (In batch) Run iTOL and save the figures in directories
