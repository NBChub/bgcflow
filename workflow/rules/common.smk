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

# set up sample
samples = pd.read_csv(config["samples"], sep="\t").set_index("strain", drop=False)
samples.index.names = ["strain_id"]

##### Wildcard constraints #####
wildcard_constraints:
    strains="|".join(samples.index),

##### Helper functions #####

STRAINS = samples.strain.to_list()
FASTA = {k: F"data/raw/internal/fasta/{v}.fasta" for (k,v) in samples.strain.to_dict().items()}