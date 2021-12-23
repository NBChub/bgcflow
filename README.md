# BGCflow
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics.svg?branch=master)](https://travis-ci.org/snakemake-workflows/snakemake-bgc-analytics)

Snakemake workflow to systematically analyze Biosynthetic Gene Clusters from a collection of genomes (pangenomes) from internal &amp; public genome datasets.

## Workflow overview
![dag](workflow/report/images/rulegraph.svg)
## Usage
### Step 1: Clone the workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis. _Make sure to have the right access / SSH Key._

    git clone git@github.com:NBChub/bgcflow.git
    cd bgcflow

### Step 2: Configure the workflow
#### 2.1 Setting Up Your Project Information
Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure your `project` and the workflow execution.

In the `config.yaml`, each `project` starts with "`-`" and must contain the name of your project (`name`) and the location of the sample file (`samples.csv`):

```yaml
projects:
  - name: example
    samples: config/samples.csv
```

Ideally, a project is a set of genomes of certain clade (pangenome) and can be annotated using a set of well-known reference genomes.

See [further configuration](#further-configuration) for more details.

#### 2.2 Setting Up Your Samples Information
The parameter `samples` denote the location of your `samples.csv` which specify the strains to analyse.

`samples.csv` example:

| genome_id       | source | organism                        | genus        | species | strain     |
|----------------:|-------:|--------------------------------:|-------------:|--------:| ----------:|
| NBC_01270.1     | custom | Streptomyces sp. NBC_01270      | Streptomyces | sp.     | NBC001270  |
| GCF_000359525.1 | ncbi   | Streptomyces albus strain J1074 | Streptomyces | albus   | J1074      |

Further formatting rules are defined in the `workflow/schemas/` folder.

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

    snakemake --use-conda --cores {number} --keep-going

Check you job DAG by executing:

    snakemake --dag | dot -Tsvg > workflow/report/images/dag.svg

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

## Further configuration
### Custom Prokka database
You can add an optional parameters: `prokka-db`, which refer to the location of a `.csv` file containing a list of your custom reference genomes for [`prokka`](https://github.com/tseemann/prokka#option---proteins) annotation:
```yaml
projects:
  - name: example
    samples: config/samples.csv
    prokka-db: config/prokka-db.csv
```

The file `prokka-db.csv` should contain a list of high quality annotated genomes that you would like to use to prioritise prokka annotations. 

`prokka-db.csv` example for Actinomycete group:

| Accession       | Strain Description             |
|----------------:|-------------------------------:|
| GCA_000203835.1 | Streptomyces coelicolor A3(2)  |
| GCA_000196835.1 | Amycolatopsis mediterranei U32 |

### Running multiple projects
You can have multiple projects running by starting a new line of project information with "`-`":

```yaml
# Project 1
  - name: example
    samples: config/samples.csv
    prokka-db: config/prokka-db.csv
# Project 2
  - name: dummy_test
    samples: config/samples_test.csv
```
Note that each `project` must have unique `name` and `samples` value.
### Choosing which analysis to run
You can choose which analysis to run by setting the parameter value to `TRUE` or `FALSE`:
```yaml
rules:
  bigscape: TRUE
  mlst: TRUE
  refseq_masher: TRUE
  seqfu: TRUE
  eggnog: FALSE
  ```
### Setting custom resources/databases folder
By default, the resources folder containing software and database dependencies are stored in the `resources/` directory. 

If you already have the resources folder somewhere else in your local machine, you can tell the workflow about their locations:

```yaml
resources_path: 
  antismash_db: $HOME/your_local_directory/antismash_db
  eggnog_db: $HOME/your_local_directory/eggnog_db
  BiG-SCAPE: $HOME/your_local_directory/BiG-SCAPE
```
