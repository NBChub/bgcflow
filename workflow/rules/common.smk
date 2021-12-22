import os
import pandas as pd
import yaml, json, sys, itertools
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("6.7.0")
__version__ = "0.1.0"


##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


##### extract project information #####
def extract_project_information():
    """
    Wrapper to extract variables from projects in config.yaml.
    Under development to accomodate multiple projects running in one snakemake run.
    """
    class ConfigError(Exception):
        """Raised when config.yaml does not satisty requirements"""
        pass

    def merge_nested_list(nested_list):
        """Merge nested list or dict values / keys"""
        return list(itertools.chain(*nested_list))
    
    # load information from config
    projects = pd.DataFrame(config["projects"])
    
    # check validity of {sample}.csv Value should be unique
    check_duplicates = projects[projects.samples.duplicated()]
    if len(check_duplicates) > 0:
        raise ConfigError(f"Project: {check_duplicates.name.to_list()} input file name: {check_duplicates.samples.to_list()} is not unique. Check your config.yaml configuration.")

    samples = []
    for i in projects.index:
        df1 = pd.read_csv(projects.loc[i, "samples"])
        df1["sample_paths"] = projects.loc[i, "samples"]
        df1["name"] = projects.loc[i, "name"]
        df1 = df1.set_index('genome_id', drop=False)
        samples.append(df1)
    df_samples = pd.concat(samples, axis=0)

    # check validity of genome_ids. Value should be unique.
    check_duplicates = df_samples[df_samples.genome_id.duplicated()]
    if len(check_duplicates) > 0:
        raise ConfigError(f"Strain ids in: {check_duplicates.sample_paths.to_list()} are not unique. Check your config.yaml configuration.")

    prokka_db = []
    for i in projects.index:
        try:
            df2 = pd.read_csv(projects.loc[i, "prokka-db"])
            df2["name"] = projects.loc[i, "name"]
            df2 = df2.set_index('Accession', drop=False)
            prokka_db.append(df2)
            df_prokka_db = pd.concat(prokka_db, axis=0).reset_index(drop=True)
        except (ValueError, KeyError) as e:
            df_prokka_db = pd.DataFrame(columns=["Accession"])
            pass
    
    return df_samples, df_prokka_db

DF_SAMPLES, DF_PROKKA_DB = extract_project_information()


##### Helper functions #####
PROJECT_IDS = DF_SAMPLES.name.unique()
STRAINS = DF_SAMPLES.genome_id.to_list()
CUSTOM = DF_SAMPLES[DF_SAMPLES.source.eq("custom")].genome_id.to_list()
NCBI = DF_SAMPLES[DF_SAMPLES.source.eq("ncbi")].genome_id.to_list()
PATRIC = DF_SAMPLES[DF_SAMPLES.source.eq("patric")].genome_id.to_list()
PROKKA_DB = DF_PROKKA_DB.Accession.to_list()
SAMPLE_PATHS = DF_SAMPLES.sample_paths.unique()


##### Helper for lambda functions #####
def get_project_only_strains(wildcards, df_samples=DF_SAMPLES):
    """
    Given a project name, extract the corresponding strain ids
    """
    PROJECT_STRAINS = {name : DF_SAMPLES[DF_SAMPLES.name.eq(name)].genome_id.to_list() for name in PROJECT_IDS}
    output = PROJECT_STRAINS[wildcards.name]
    return output

# prokka #
def get_prokka_db_accessions(wildcards, df_prokka_db=DF_PROKKA_DB):
    """
    Given a project name, find which accessions to create a prokka reference gbff
    """
    accession = df_prokka_db[df_prokka_db["name"] == wildcards.name].Accession.to_list()
    output = [f"resources/prokka_db/gbk/{acc}.gbff" for acc in accession]
    return output

def get_prokka_refdb(wildcards, df_samples=DF_SAMPLES):
    """
    Given a strain id, find which prokka db to use
    """
    name = df_samples[df_samples["genome_id"] == wildcards.strains].name.values
    if name in DF_PROKKA_DB.name.unique():
        output = f"--proteins resources/prokka_db/reference_{name[0]}.gbff"
    else:
        output = ""
    return output

def get_prokka_sample_path(wildcards, df_samples=DF_SAMPLES):
    """
    Given a project name, find the corresponding sample file
    """
    output = df_samples[df_samples["genome_id"] == wildcards.name].samples_path.values
    return output

# bigscape #
def get_bigscape_inputs(name, version, df_samples=DF_SAMPLES):
    """
    Given a project name, find the corresponding sample file
    """
    selection = df_samples[df_samples["name"] == name].genome_id.values
    output = [f"data/interim/antismash/{version}/{s}/{s}.gbk" for s in selection]
    return output


##### Wildcard constraints #####
wildcard_constraints:
    strains="|".join(STRAINS),
    ncbi="|".join(NCBI),
    custom="|".join(CUSTOM),
    patric="|".join(PATRIC),
    name="|".join(PROJECT_IDS),


##### Get dependency versions #####
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
                "roary" : expand("data/interim/roary/{name}", name=PROJECT_IDS),
                "bigscape" : expand("data/interim/bigscape/{name}_antismash_{version}/index.html", version=dependency_version["antismash"], name=PROJECT_IDS),
                "seqfu" : "data/processed/tables/df_seqfu_stats.csv"
                }
    
    # get keys from config
    opt_rules = config["rules"].keys()

    # if values are TRUE add output files to rule all
    final_output = [rule_dict[r] for r in opt_rules if config["rules"][r]]
     
    return final_output


##### Custom Resource Directory #####
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