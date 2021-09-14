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

## watermark snakemake version to samples
for i in samples.index:
    samples.loc[i, "strain"] = str(samples.loc[i, "strain"]) + __version__

validate(samples, schema="../schemas/samples.schema.yaml")

# set up units
units = pd.read_table(config["units"], dtype=str).set_index(
    ["strain"], drop=False
)

## watermark snakemake version to units
for i in units.index:
    units.loc[i, "strain"] = str(units.loc[i, "strain"]) + __version__

validate(units, schema="../schemas/units.schema.yaml")

#units.index = units.index.set_levels(
#    [i.astype(str) for i in units.index.levels]
#)  # enforce str in index

##### Wildcard constraints #####
wildcard_constraints:
#    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),

##### Helper functions #####


STRAINS = samples.strain.to_list()
ILLUMINA = {k: v for (k,v) in units.illumina_reads.to_dict().items()}
NANOPORE = {k: v for (k,v) in units.nanopore_reads.to_dict().items()}
