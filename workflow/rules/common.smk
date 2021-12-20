import os
import pandas as pd
import yaml, json
import sys
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

##### dependency versions #####
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

##### Customizable Analysis #####
def get_final_output():
    """
    Generate final output for rule all given a TRUE value in config["rules"]
    """
    # dictionary of rules and its output files
    rule_dict = {"mlst" : expand("data/interim/mlst/{strains}_ST.csv", strains = STRAINS),
                "eggnog" : expand("data/interim/eggnog/{strains}/", strains = STRAINS),
                "refseq_masher" : expand("data/interim/refseq_masher/{strains}_masher.csv", strains = STRAINS),
                "automlst_wrapper" : "data/interim/automlst_wrapper/raxmlpart.txt.treefile",
                "roary" : "data/interim/roary/all",
                "bigscape" : expand("data/interim/bigscape/antismash_{version}/index.html", version=dependency_version["antismash"]),
                "seqfu" : "data/processed/tables/df_seqfu_stats.csv"
                }
    
    # get keys from config
    opt_rules = config["rules"].keys()

    # if values are true add output files to rule all
    final_output = [rule_dict[r] for r in opt_rules if config["rules"][r]]
    return final_output

##### Custome Resource Directory #####
def custom_resource_dir():
    """
    Generate symlink for user defined resources location
    """
    resource_dbs = config["resources_path"]
    sys.stderr.write(f"Checking for user-defined local resources...\n")
    for r in resource_dbs.keys():
        # check for default path
        path = Path(resource_dbs[r])
        if resource_dbs[r] == f"resources/{r}":
            pass 
        # check for user-defined external resources
        elif path.exists():
            try:    
                slink = Path(f"resources/{r}")
                existing_path = Path.readlink(slink)
                # check if symlink for extrenal path is already generated
                if existing_path == path:
                    sys.stderr.write(f"- Symlink for {r} already exists at: {existing_path}\n")
                # update symlink because new path is given
                else:
                    slink.unlink()
                    slink.symlink_to( path )
                    updated_path = Path.readlink(Path(slink)) 
                    sys.stderr.write(f"- Updating symlink for {r} from: {existing_path} to: {updated_path}\n")
            # generate a new symlink
            except FileNotFoundError:
                sys.stderr.write(f"- Generating symlink for {r} from: {path}\n")                    
                slink.symlink_to( path )
        # raise an Error if external path not found
        else:
            raise FileNotFoundError(f"Error: User-defined resource {r} at {path} does not exist. Check the config.yaml and provide the right path for resource {r} or change it to the default path: resources/{r}\n")
    return 
custom_resource_dir()