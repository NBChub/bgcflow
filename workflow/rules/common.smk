import os
import pandas as pd
import yaml, json
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("6.7.0")
__version__ = "0.1.0"

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

# dependency versions
def get_dependency_version(dep, dep_key):
    """
    return dependency version tags given a dictionary (dep) and its key (dep_key)
    """
    with open(dep[dep_key]) as file:
        result = []
        documents = yaml.full_load(file)
        for i in documents["dependencies"]:
            if i.startswith(dep_key):
                result = i.split("=")[-1]
    return str(result)

def write_dependecies_to_json(dep, outfile):
    """
    write dependency version to a json file
    """
    with open(outfile, "w") as file:
        dv = {}
        for ky in dependencies.keys():
            vr = get_dependency_version(dependencies, ky)
            dv[ky] = vr
        json.dump(dv, file, indent=2,)
        file.close()
    return dv

# list of the main dependecies used in the workflow
dependencies = {"antismash" : r"workflow/envs/antismash.yaml",
                "prokka": r"workflow/envs/prokka.yaml",
                "mlst" : r"workflow/envs/mlst.yaml",
                "eggnog-mapper" : r"workflow/envs/eggnog.yaml",
                "roary" : r"workflow/envs/roary.yaml",
                "refseq_masher" : r"workflow/envs/refseq_masher.yaml",
                "seqfu" : r"workflow/envs/seqfu.yaml"
                }

dependency_version = write_dependecies_to_json(dependencies, "workflow/report/dependency_versions.json")