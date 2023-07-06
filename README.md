# BGCFlow
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.14.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![PEP compatible](https://pepkit.github.io/img/PEP-compatible-green.svg)](https://pep.databio.org)

BGCFlow is a systematic workflow for the analysis of biosynthetic gene clusters across large collections of genomes (pangenomes) from internal &amp; public datasets.

## Publication
> Matin Nuhamunada, Omkar S. Mohite, Patrick V. Phaneuf, Bernhard O. Palsson, and Tilmann Weber. (2023). BGCFlow: Systematic pangenome workflow for the analysis of biosynthetic gene clusters across large genomic datasets. bioRxiv 2023.06.14.545018; doi: [https://doi.org/10.1101/2023.06.14.545018](https://doi.org/10.1101/2023.06.14.545018)

## Quick Start
A quick and easy way to use BGCFlow using [`bgcflow_wrapper`](https://github.com/NBChub/bgcflow_wrapper).

1. Create a conda environment and install the [BGCFlow python wrapper](https://github.com/NBChub/bgcflow_wrapper) :

```bash
# create and activate a new conda environment
conda create -n bgcflow pip -y
conda activate bgcflow

# install BGCFlow wrapper
pip install git+https://github.com/NBChub/bgcflow_wrapper.git
```

2. Deploy and run BGCFlow:
```bash
# Deploy and run BGCFlow
BGCFLOW_PATH="<change this to your desired path for BGCFlow>"
bgcflow clone $BGCFLOW_PATH # clone BGCFlow to BGCFLOW_PATH
cd $BGCFLOW_PATH # move to BGCFLOW_PATH
bgcflow init # initiate BGCFlow config and examples from template
bgcflow run -n # do a dry run, remove the flag "-n" to run the example dataset
```

See [`README.md`](https://github.com/NBChub/bgcflow_wrapper) for more details about [`bgcflow_wrapper`](https://github.com/NBChub/bgcflow_wrapper).

[![asciicast](https://asciinema.org/a/595149.svg)](https://asciinema.org/a/595149)

## Workflow overview
The main Snakefile workflow comprises various pipelines for data selection, functional annotation, phylogenetic analysis, genome mining, and comparative genomics for Prokaryotic datasets.

![dag](workflow/report/images/rulegraph_annotated.png)

Available pipelines in the main Snakefile can be checked using the following command:
```
bgcflow pipelines
```

## Usage
### Step 1: Clone the workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis. _Make sure to have the right access / SSH Key._

    git clone git@github.com:NBChub/bgcflow.git
    cd bgcflow

> **TIPS** - **`bgcflow_wrapper`** equivalent:
>
> ```bash
> bgcflow clone bgcflow
> ```

### Step 2: Configure the workflow
Configure the workflow according to your needs by editing the files in the `config/` folder.

#### 2.1 Using template example
An example of the configuration files is provided in the `.examples` folder.

If you have a fresh copy of BGCFlow, you can initiate config and examples using by copying the necessary files to `config/` folder:
```shell
cp .examples/_config_example.yaml config/config.yaml
```

The above command will create a new file in `config/config.yaml`. You can adjust the `config.yaml` to configure your `project` and the workflow execution.

> **TIPS** - **`bgcflow_wrapper`** equivalent
>
> ```bash
> bgcflow init --bgcflow_dir bgcflow
> ```

#### 2.2 Configure your project
##### 2.2.1 PEP Format
As of BGCFlow version `>=0.4.0`, projects are now configured as a [Portable Encapsulated Project (PEP)](http://pep.databio.org/en/latest/). In the main `config/config.yaml`, each `project` starts with "`-`" and the variable `name` which points to a PEP config file.

```yaml
projects:
  - name: .examples/_pep_example/project_config.yaml
```
See [project_config.yaml](.examples/_pep_example/project_config.yaml) for an example of a PEP formatted project.

> **TIPS** - Initiate a project using **`bgcflow_wrapper`**:
>
> ```bash
> bgcflow init --project MY_PROJECT --bgcflow_dir bgcflow
> ```

##### 2.2.1 BGCFlow Format
A project can also be configured as previously described in BGCFlow version `<=0.3.3`. In the main `config/config.yaml`, each `project` starts with "`-`" and must contain the name of your project (`name`), the location of the sample file (`samples.csv`), and a rule configuration file (`project_config.csv`):

```yaml
projects:
  - name: example
    samples: .examples/_genome_project_example/samples.csv
    rules: .examples/_genome_project_example/project_config.yaml
```
Note that the location of the sample file and the rule configuration file is relative to your `bgcflow` directory.

Ideally, you can organize a project as a set of genomes from a certain clade (pangenome).

See [further configuration](#further-configuration) for more details.

#### 2.2 Setting Up Your Samples Information
The variable `sample_table` (PEP) or `samples` denote the location of your `.csv` file which specifies the genomes to analyze. Note that you can name the file anything as long as you define it in the `config.yaml`.

Example: `samples.csv`

| genome_id       | source | organism                        | genus        | species | strain     |closest_placement_reference|
|----------------:|-------:|--------------------------------:|-------------:|--------:| ----------:|--------------------------:|
| GCF_000359525.1 | ncbi   |                                 |              |         | J1074      |                           |
| 1223307.4       | patric | Streptomyces sp. PVA 94-07      | Streptomyces | sp.     | PVA 94-07  | GCF_000495755.1           |
| P8-2B-3.1       | custom | Streptomyces sp. P8-2B-3        | Streptomyces | sp.     | P8-2B-3    |                           |

Columns description:
- **`genome_id`** _[required]_:  The genome accession ids (with genome version for `ncbi` and `patric` genomes). For `custom` fasta file provided by users, it should refer to the fasta file names stored in the `data/raw/fasta/` directory with `.fna` extension. **Example:** genome id P8-2B-3.1 refers to the file `data/raw/fasta/P8-2B-3.1.fna`.
- **`source`** _[required]_: Source of the genome to be analyzed choose one of the following: `custom`, `ncbi`, `patric`. Where:
  - `custom`: for user-provided genomes (`.fna`) in the `data/raw/fasta` directory with genome ids as filenames
  - `ncbi`: for list of public genome accession IDs that will be downloaded from the NCBI refseq (GCF...) or genbank (GCA...) database
  - `patric`: for list of public genome accession IDs that will be downloaded from the PATRIC database
- `organism` _[optional]_: name of the organism that is the same as in the fasta header
- `genus` _[optional]_ : genus of the organism. Ideally identified with GTDBtk.
- `species` _[optional]_: species epithet (the second word in a species name) of the organism. Ideally identified with GTDBtk.
- `strain` _[optional]_ : strain id of the organism
- `closest_placement_reference` _[optional]_: if known, the closest NCBI genome to the organism. Ideally identified with GTDBtk.

Further formatting rules are defined in the `workflow/schemas/` folder.

#### 2.3 Rules: Choosing which analysis to run
In each projects, you can choose which analysis to run by setting the parameter value in the [`project_config.yaml`](.examples/_genome_project_example/project_config.yaml) to `TRUE` or `FALSE`:
```yaml
rules:
  bigscape: TRUE
  mlst: TRUE
  refseq_masher: TRUE
  seqfu: TRUE
  eggnog: FALSE
  ```
> **TIPS** - Finding available rules with **`bgcflow_wrapper`**
>
> ```bash
> bgcflow pipelines --bgcflow_dir bgcflow
> ```

> **TIPS** - Find out rule description with **`bgcflow_wrapper`**
>
> ```bash
> bgcflow pipelines --describe bigscape --bgcflow_dir bgcflow
> ```

See [List of Configurable Features](##List-of-Configurable-Features) for more details.
### Step 3: Install Snakemake & BGCFlow environment

Installing Snakemake using [Mamba](https://github.com/mamba-org/mamba) is advised. In case you don’t use [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) you can always install [Mamba](https://github.com/mamba-org/mamba) into any other Conda-based Python distribution with:

    conda install -n base -c conda-forge mamba

You can use [`bgcflow_wrapper`](https://github.com/NBChub/bgcflow_wrapper) environment from [Quick Start](#Quick-Start) or install BGCFlow environment which contain Snakemake (`version 7.14.0`) and other dependencies with:

```bash
# create and activate a new conda environment
conda create -n bgcflow pip -y
conda activate bgcflow

# install BGCFlow wrapper
pip install git+https://github.com/NBChub/bgcflow_wrapper.git
```

See the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for installation details.

### Step 4: Execute workflow

Activate the conda environment:

    conda activate bgcflow

Test your configuration by performing a dry-run via:

    snakemake --use-conda -n

Execute the workflow locally via:

    snakemake --use-conda --cores {number} --keep-going

Check you job DAG by executing:

    snakemake --dag | dot -Tsvg > workflow/report/images/dag.svg

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

## Further configuration
### Custom Prokka database
You can add an optional parameter: `prokka-db`, which refers to the location of a `.csv` file containing a list of your custom reference genomes for [`prokka`](https://github.com/tseemann/prokka#option---proteins) annotation:
```yaml
projects:
  - name: example
    samples: config/samples.csv
    prokka-db: config/prokka-db.csv
```

The file `prokka-db.csv` should contain a list of high-quality annotated genomes that you would like to use to prioritize prokka annotations.

`prokka-db.csv` example for Actinomycete group:

| Accession       | Strain Description             |
|----------------:|-------------------------------:|
| GCA_000203835.1 | Streptomyces coelicolor A3(2)  |
| GCA_000196835.1 | Amycolatopsis mediterranei U32 |

### Taxonomic Placement
The workflow will prioritize user-provided taxonomic placement by adding an optional parameter: `gtdb-tax`, which refers to a similar GTDB-tk summary file, but only the "user_genome" and "classification" columns are required.

`gtdbtk.bac120.summary.tsv` example:

| user_genome | classification                                                                                                                           |
|------------:|---------------------------------------------------------------------------------------------------------------------------------------:|
| P8-2B-3.1   | d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Streptomycetales;f__Streptomycetaceae;g__Streptomyces;s__Streptomyces albidoflavus |

If these are not provided, the workflow will use the `closest_placement_reference` columns in the sample file (see above). Note that the value must be a valid genome accession in the latest GTDB release (currently R202), otherwise, it will raise an error.

If this information is not provided, then the workflow will guess the taxonomic placement by:
1. If the `source` is `ncbi`, it will try to find the accession via GTDB API. If it doesn't find any information then,
2. It will use the `genus` table and find the parent taxonomy via GTDB API, which then results in `_genus_ sp.` preceded by the matching parent taxonomy.
3. If both option does not find any taxonomic information, then it will return empty taxonomic values.

### Running multiple projects
You can have multiple projects running by starting a new line of project information with "`-`":

```yaml
projects:
# Project 1
  - name: example
    samples: config/samples.csv
    prokka-db: config/prokka-db.csv
# Project 2
  - name: example_2
    samples: config/samples_2.csv
```
Note that each `project` must have a unique `name` and `samples` value.

### Setting custom resources/databases folder
By default, the resources folder containing software and database dependencies is stored in the `resources/` directory.

If you already have the resources folder somewhere else in your local machine, you can tell the workflow about their locations:

```yaml
resources_path:
  antismash_db: $HOME/your_local_directory/antismash_db
  eggnog_db: $HOME/your_local_directory/eggnog_db
  BiG-SCAPE: $HOME/your_local_directory/BiG-SCAPE
```
## List of Configurable Features
Here you can find rules keywords that you can run within BGCflow.
| Keywords | Description | Links |
|:---------| :------------- | :------------------------- |
| seqfu | Returns contig statistics of the genomes | [SeqFu](https://github.com/telatin/seqfu2)|
| mlst | Returns genome classification within multi-locus sequence types (STs) | [mlst](https://github.com/tseemann/mlst) |
| refseq_masher | Identify the closest 10 NCBI Refseq genomes | [RefSeq Masher](https://github.com/phac-nml/refseq_masher) |
| mash | Calculate genomic distance using MAST | [MASH](https://github.com/marbl/Mash) |
| fastani | Calculate nucleotide distance using fastANI | [fastANI](https://github.com/ParBLiSS/FastANI) |
| checkm | Assess genome quality | [CheckM](https://github.com/Ecogenomics/CheckM) |
| gtdbtk | Identify taxonomy of genomes using GTDB-toolkit | [GTDBTk](https://github.com/Ecogenomics/GTDBTk) |
| prokka-gbk | Returns annotated `.gbk` files | [Prokka](https://github.com/tseemann/prokka) |
| diamond | Create diamond database for alignment | [DIAMOND](https://github.com/bbuchfink/diamond) |
| antismash-summary | Summary of BGCs statistics | [antiSMASH](https://github.com/antismash/antismash) |
| antismash-zip | Returns zipped antiSMASH result | [antiSMASH](https://github.com/antismash/antismash) |
| query_bigslice | Query BGCs with BiG-FAM db* | [BiG-SLICE](https://github.com/medema-group/bigslice) |
| bigscape | Build Gene Cluster Families with BiG-SCAPE | [BiG-SCAPE](https://git.wageningenur.nl/medema-group/BiG-SCAPE) |
| bigslice | Build Gene Cluster Families with BiG-SLICE | [BiG-SLICE](https://github.com/medema-group/bigslice) |
| automlst_wrapper | Build a phylogenomic tree with autoMLST wrapper | [autoMLST-wrapper](https://github.com/KatSteinke/automlst-simplified-wrapper), [autoMLST](https://bitbucket.org/ziemertlab/automlst/src/master/) |
| roary | Build Pangenome | [Roary](https://github.com/sanger-pathogens/Roary) |
| eggnog | Functional annotation with EggNOG-mapper | [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) |
| deeptfactor | Prediction of transcription factors with DeepTFactor | [DeepTFactor](https://bitbucket.org/kaistsystemsbiology/deeptfactor) |
| roary++ | Apply multiple tools together with Roary pangenome (diamond, automlst_wrapper, eggnog, deeptfactor) | [Roary](https://github.com/sanger-pathogens/Roary)  |
| cblaster-genome | Generate cblaster databases for genomes in project | [cblaster](https://github.com/gamcil/cblaster)  |
| cblaster-bgcs | Generate cblaster databases for bgcs in project | [cblaster](https://github.com/gamcil/cblaster)  |

## Using snakemake profiles for further configurations
When using different machines, you can, for example, adapt the number of threads required for each rule using a Snakemake profile. An example is given in [`config/examples/_profile_example/config.yaml`](config/examples/_profile_example/config.yaml):
```yaml
set-threads:
  - antismash=4
  - arts=4
  - bigscape=32
  - bigslice=16
```

You can use run a snakemake job with the above profile with:
```bash
snakemake --profile config/examples/_profile_example/ --use-conda -c $N -n # remove the dry-run parameters "-n" for the actual run
```
Or also with a defined `config` file:
```bash
snakemake --configfile config/examples/_config_example.yaml --profile config/examples/_profile_example/ --use-conda -c $N -n # remove the dry-run parameters "-n" for the actual run
```
