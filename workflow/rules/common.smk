import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# set up sample for default case with fasta files provided
df_samples = pd.read_csv(config["samples"]).set_index("genome_id", drop=False)
df_samples.index.names = ["genome_id"]

# set up custom prokka database
df_prokka_db = pd.read_csv(config["prokka-db"]).set_index("Accession", drop=False)
df_prokka_db.index.names = ["Accession"]

##### Helper functions #####
STRAINS = df_samples.genome_id.to_list()
CUSTOM = df_samples[df_samples.source.eq("custom")].genome_id.to_list()
NCBI = df_samples[df_samples.source.eq("ncbi")].genome_id.to_list()
PATRIC = df_samples[df_samples.source.eq("patric")].genome_id.to_list()
PROKKA_DB = df_prokka_db.Accession.to_list()

# Helper for lamda function
SAMPLES = {idx : idx for idx in STRAINS}

##### Wildcard constraints #####
wildcard_constraints:
    strains="|".join(STRAINS),
    ncbi="|".join(NCBI),
    custom="|".join(CUSTOM),
    patric="|".join(PATRIC),
    prokka_db="|".join(PROKKA_DB)