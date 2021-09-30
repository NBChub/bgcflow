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
df_samples = pd.read_csv(config["samples"], sep="\t").set_index("genome_id", drop=False)
df_samples.index.names = ["genome_id"]

##### Wildcard constraints #####
wildcard_constraints:
    strains="|".join(df_samples.index),

##### Helper functions #####
STRAINS = df_samples.genome_id.to_list()
FASTA = dict()

fasta_input_dir = "data/raw/fasta/"
for genome_id in STRAINS:
    if df_samples.loc[genome_id, 'source'] == 'custom':
        if genome_id + '.fna' in os.listdir():
            FASTA[genome_id] = os.path.join(fasta_input_dir, genome_id + '.fna')
        # else:
        #     logging('Error: Genome does not have corresponding fasta file in raw data folder:', genome_id)
        
    # elif df_samples.loc[genome_id, 'source'] == 'ncbi':
    #    then activate NCBI rule and download fna in fasta directory
    #           once the fna is downloaded in right folder 
    #           then FASTA[genome_id] = os.path.join(fasta_input_dir, genome_id + '.fna')
    # elif source has azure id
    #    then activate az copy and pull the fna from azure
    #           once the fna is downloaded in right folder 
    #           then FASTA[genome_id] = os.path.join(fasta_input_dir, genome_id + '.fna')
    # else
    #    undefined source found for the genome use one of the above sources
