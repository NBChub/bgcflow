import os
import numpy as np
import pandas as pd
import yaml, json, sys, itertools, hashlib
from snakemake.utils import validate
from snakemake.utils import min_version
from pathlib import Path

min_version("7.14.0")
__version__ = "0.4.0"

container: "docker://continuumio/miniconda3:4.12.0"

##### TABLE OF CONTENTS #####
# This .smk file helps the rules interact with the user config and contains 
# many helper scripts and python functions. Execution of the functions should
# be done in the Snakefile whenever possible. This will help in creating
# custom sub-workflows (Snakefile) in the future.
# In general, the structure of a workflow can be divided into these steps:
#   1. Load config file
#   2. Extract project information (called via Snakefile)
#   3. Generate wildcard constants (via Snakefile)
#   4. Constraints wildcards (via Snakefile)
#   5. Helper lambda functions for rules I/O (called by specific rules (.smk))
#   6. Get dependency versions
#   7. Customize final output based on config["rule"] values (called via Snakefile)
#   8. Set up custom resource directory provided in config["resources_path"] (called via Snakefile)


##### 1. Load config #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

sys.stderr.write(f"This is BGCflow version {__version__}.\n\n")

##### 2. Extract project information #####
# The function extract_project_information() returns objects necessary for steps 3 and 4
# The function is called in Snakefile where wildcards are extracted and defined

def extract_project_information(config):
    """
    Wrapper to extract variables from projects in config.yaml.
    """

    class ConfigError(Exception):
        """Raised when config.yaml does not satisty requirements"""
        pass


    def hash_prokka_db(prokka_db_path):
        """
        Read a csv containing a list of accession ids and generate a hash identifier.
        The list must have "Accession" as the column name.
        """
        df = pd.read_csv(Path(prokka_db_path))
        assert "Accession" in df.columns
        column = df.loc[:, "Accession"].sort_values().reset_index(drop=True)
        hash_value = hashlib.sha1(pd.util.hash_pandas_object(column).values).hexdigest()

        # Dictionary
        hash_object = {hash_value : column.to_list()}

        # File mapping
        file_map = {str(prokka_db_path) : hash_value}
        return hash_object, file_map


    def create_prokka_db_tables(prokka_db_path):
        """
        Build a json dictionaries of prokka references and files from config information
        """
        prokka_db_table = {}
        prokka_db_map = {}

        if Path(prokka_db_path).suffix == ".csv":
            hash_value, hash_map = hash_prokka_db(prokka_db_path)
            prokka_db_table.update(hash_value)
            prokka_db_map.update(hash_map)
        elif Path(prokka_db_path).suffix == ".gbff":
            prokka_db_map.update(prokka_db_path)
        else:
            raise ValueError('Prokka-db should be a csv table containing valid genome accession ids or a collection of annotated genomes in .gbff format')
        return prokka_db_table, prokka_db_map


    def find_conflicting_samples(df_samples):
        """
        Prune identical genome_ids across projects
        """
        filtered_columns = ['sample_paths',
                            'prokka-db',
                            'gtdb_paths',
                            'name']

        df_filtered = df_samples.copy().drop(columns=filtered_columns)
        dropped_ids = df_filtered[df_filtered.duplicated()].index
        df_filtered = df_filtered.drop_duplicates()
        if df_filtered.index.has_duplicates:
            duplicates = df_filtered[df_filtered.index.duplicated()].index
            duplicates = df_samples.loc[duplicates]
            for genome_id in duplicates['genome_id'].unique():
                print(f"WARNING: Samples {genome_id} in projects {duplicates.loc[genome_id, 'name'].to_list()} is not identical.", file=sys.stderr)
            print("Each samples with the same genome_id should have identical annotation and taxonomy. Please use different ids.", file=sys.stderr)
            print(duplicates, file=sys.stderr)
            raise ConfigError

        for ids in df_filtered.index:
            for col in filtered_columns:
                value = df_samples.loc[ids, col]
                if type(value) == str:
                    df_filtered.at[ids, col] = [value]
                else:
                    try:
                        df_filtered.at[ids, col] = value.unique()
                    except AttributeError as e:
                        pass

        return df_filtered

    # load information from config
    print(f"Step 1. Extracting project information from config...\n", file=sys.stderr)
    projects = config["projects"]
    df_projects = pd.DataFrame(projects).set_index('name', drop=False)

    df_samples = []
    prokka_db_table = {}
    prokka_db_map = {}

    for num, p in enumerate(projects):
        # load samples
        print(f"Step 2.{num+1} Getting sample information for project: {p['name']}", file=sys.stderr)
        df_sample = pd.read_csv(p['samples']).set_index("genome_id", drop=False)
        df_sample.loc[:, 'name'] = p['name']
        df_sample.loc[:, 'sample_paths'] = p['samples']
        df_sample.loc[:, 'prokka-db'] = np.nan
        df_sample.loc[:, 'gtdb_paths'] = np.nan

        df_sample.loc[:, 'classification'] = np.nan
        df_sample.loc[:, 'classification_source'] = np.nan
        # load prokka_db
        if 'prokka-db' in p.keys():
            print(f" - Found user-provided reference genomes for Prokka annotation", file=sys.stderr)
            prokka_db_path = p['prokka-db']
            df_sample.loc[:, 'prokka-db'] = prokka_db_path
            prokka_table, prokka_map = create_prokka_db_tables(prokka_db_path)
            prokka_db_table.update(prokka_table)
            prokka_db_map.update(prokka_map)

            df_sample["ref_annotation"] = list(prokka_table.keys())[0]
        if 'gtdb-tax' in p.keys():
            print(f" - Found user-provided taxonomic information\n", file=sys.stderr)
            df_sample.loc[:, 'gtdb_paths'] = p['gtdb-tax']
            df_gtdb = pd.read_csv(p['gtdb-tax'], sep="\t").set_index("user_genome")
        else:
            df_gtdb = pd.DataFrame()

        for i in df_sample.index:
            if i in df_gtdb.index:
                df_sample.loc[i, "classification"] = df_gtdb.loc[i, "classification"]
                df_sample.loc[i, "classification_source"] = "user_provided"
            elif not pd.isnull(df_sample.loc[i, "closest_placement_reference"]):
                df_sample.loc[i, "classification"] = df_sample.loc[i, "closest_placement_reference"]
                df_sample.loc[i, "classification_source"] = "ncbi"
            elif not pd.isnull(df_sample.loc[i, "genus"]):
                df_sample.loc[i, "classification"] = df_sample.loc[i, "closest_placement_reference"]
                df_sample.loc[i, "classification_source"] = "user_provided_genus"
            elif df_sample.loc[i, "source"] in ['ncbi', 'patric']:
                df_sample.loc[i, "classification_source"] = df_sample.loc[i, "source"]
        df_samples.append(df_sample)

    df_samples = pd.concat(df_samples)

    print(f"Step 3. Merging genome_ids across projects...\n", file=sys.stderr)
    df_samples = find_conflicting_samples(df_samples)

    print(f"Step 4. Checking validity of samples using schemas..\n", file=sys.stderr)
    filtered_columns = ['sample_paths',
                        'prokka-db',
                        'gtdb_paths',
                        'name']
    validate(df_samples.drop(columns=filtered_columns).fillna(""), schema="schemas/samples.schema.yaml")

    return df_projects, df_samples, prokka_db_table, prokka_db_map

##### 3. Generate wildcard constants
# This step is done via Snakefile. Refer to comments on Step 2 for explanation.

##### 4. Wildcard constraints
# This step is done via Snakefile. Refer to comments on Step 2 for explanation.

##### 5. Helper lambda functions for calling rules I/O #####

# seqfu.smk #
def get_fasta_inputs(name, df_samples):
    """
    Given a project name, list all corresponding strains (genome_id) fasta file
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/fasta/{s}.fna" for s in selection]
    return output

# prokka.smk #
def get_prokka_refdb(genome_id, params, df_samples, mapping_file, config=config):
    """
    Given a genome id, find which prokka-db input to use.
    params:
        - "table" - will return the corresponding prokka-db table to use
        - "file" - will return the corresponding reference gbks
        - "params" - will return prokka protein params and the corresponding file
    """

    prokka_db = df_samples.loc[genome_id, "prokka-db"][0]
    name = df_samples.loc[genome_id, "name"][0]

    if not os.path.isfile(str(prokka_db)):
        if params == "file":
            output = []
        else:
            output = ""
    elif params == "table":
        output = prokka_db
    elif params == "file":
        output = f"resources/prokka_db/{mapping_file[prokka_db]}.gbff"
    elif params == "params":
        output = f"--proteins resources/prokka_db/{mapping_file[prokka_db]}.gbff"
    else:
        sys.stderr.write(f"Second argument should be: table, file, or params.\n")
        raise
    return output

# bigscape.smk, bigslice.smk, and bgc_analytics.smk #
def get_antismash_inputs(name, version, df_samples):
    """
    Given a project name, find the corresponding sample file to use
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/antismash/{version}/{s}/{s}.gbk" for s in selection]
    return output

# roary.smk #
def get_prokka_outputs(name, df_samples, ext="gff"):
    """
    Given a project name, find the corresponding sample file to use
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/prokka/{s}/{s}.{ext}" for s in selection]
    return output

# automlst_wrapper.smk #
def get_automlst_inputs(name, df_samples):
    """
    Given a project name, find the corresponding sample file to use
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/automlst_wrapper/{name}/{s}.gbk" for s in selection]
    return output

# gtdb.smk #
def get_json_inputs(name, df_samples):
    """
    Given a project name, find the corresponding sample file to use
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    output = [f"data/interim/gtdb/{s}.json" for s in selection]
    return output

# ncbi.smk #
def get_ncbi_assembly_inputs(name, df_samples):
    """
    Given a project name, find the corresponding sample file to use
    """
    selection = [i for i in df_samples.index if name in df_samples.loc[i, "name"]]
    selection_ncbi = df_samples[df_samples["source"] == "ncbi"].genome_id.values
    output = [f"data/interim/assembly_report/{s}.json" for s in selection_ncbi]
    return output

##### 6. Get dependency versions #####
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

def get_dependencies(dep):
    """
    get dependency version
    """
    dv = {}
    for ky in dep.keys():
        vr = get_dependency_version(dep, ky)
        dv[ky] = vr
    return dv


# list of the main dependecies used in the workflow
dependencies = {"antismash" : r"workflow/envs/antismash.yaml",
                "prokka": r"workflow/envs/prokka.yaml",
                "mlst" : r"workflow/envs/mlst.yaml",
                "eggnog-mapper" : r"workflow/envs/eggnog.yaml",
                "roary" : r"workflow/envs/roary.yaml",
                "refseq_masher" : r"workflow/envs/refseq_masher.yaml",
                "seqfu" : r"workflow/envs/seqfu.yaml",
                "checkm" : r"workflow/envs/checkm.yaml",
                }

dependency_version = get_dependencies(dependencies)

##### 7. Customize final output based on config["rule"] values #####
def get_project_outputs(config, PROJECT_IDS, df_samples):
    """
    Generate outputs of a project given a TRUE value in config["rules"]
    """

    # read rule configs
    if type(config) == dict:
        pass
    else:
        with open(config, 'r') as file:
            config = yaml.safe_load(file)

    selection = [i for i in df_samples.index if PROJECT_IDS in df_samples.loc[i, "name"]]

    # dictionary of rules and its output files
    rule_dict = {"mlst" : expand("data/interim/mlst/{strains}_ST.csv", strains = selection),
                "eggnog" : expand("data/interim/eggnog/{strains}/", strains = selection),
                "refseq-masher" : expand("data/interim/refseq_masher/{strains}_masher.csv", \
                                         strains = selection),
                "mash" : expand("data/interim/mash/{name}/triangle_distance_matrix.tsv", \
                                         name=PROJECT_IDS),
                "fastani" : expand("data/interim/fastani/{name}/fastani_out.tsv", \
                                         name=PROJECT_IDS),
                "automlst-wrapper" : expand("data/processed/{name}/automlst_wrapper/final.newick", \
                                            name=PROJECT_IDS),
                "roary" : expand("data/processed/{name}/roary/df_gene_presence_binary.csv", name=PROJECT_IDS),
                "eggnog-roary" : expand("data/interim/eggnog_roary/{name}/", name=PROJECT_IDS),
                "seqfu" : expand("data/processed/{name}/tables/df_seqfu_stats.csv", name=PROJECT_IDS),
                "rnammer": "resources/rnammer_test.txt",
                "bigslice": expand("data/processed/{name}/bigslice/cluster_as_{version}/", \
                                    name = PROJECT_IDS, version=dependency_version["antismash"]),
                "query-bigslice":  expand("data/processed/{name}/bigslice/query_as_{version}/", \
                                name = PROJECT_IDS, version=dependency_version["antismash"]),
                "checkm" : expand("data/processed/{name}/tables/df_checkm_stats.csv", name = PROJECT_IDS),
                "gtdbtk" : expand("data/processed/{name}/tables/gtdbtk.bac120.summary.tsv", name = PROJECT_IDS),
                "prokka-gbk" : expand("data/processed/{name}/genbank/{strains}.gbk", strains=selection,
                                        name = PROJECT_IDS),
                "antismash-summary": expand("data/processed/{name}/tables/df_antismash_{version}_summary.csv", \
                                            name = PROJECT_IDS, version=dependency_version["antismash"]),
                "antismash-zip": expand("data/processed/{name}/antismash/{version}/{strains}.zip",
                                        name=PROJECT_IDS, strains=selection, version=dependency_version["antismash"]),
                "arts": expand("data/interim/arts/antismash-{version}/{strains}/", \
                                version=dependency_version["antismash"], strains = selection),
                "bigscape" : expand("data/processed/{name}/bigscape/for_cytoscape_antismash_{version}", \
                                     name = PROJECT_IDS, version=dependency_version["antismash"]),
                "diamond" : expand("data/processed/{name}/diamond/{name}.dmnd", name = PROJECT_IDS),
                "diamond-roary" : expand("data/processed/{name}/diamond/{name}_pangenome.dmnd", name=PROJECT_IDS),
                "deeptfactor" : expand("data/interim/deeptfactor/{strains}/", strains = selection),
                "deeptfactor-roary" : expand("data/interim/deeptfactor_roary/{name}/", name=PROJECT_IDS),
                "cblaster-genome":  expand("data/processed/{name}/cblaster/genomes/", name = PROJECT_IDS),
                "cblaster-bgc":  expand("data/processed/{name}/cblaster/bgcs/{version}/", \
                                name = PROJECT_IDS, version=dependency_version["antismash"]),
                }
    
    # get keys from config
    opt_rules = config['rules'].keys()

    # if values are TRUE add output files to rule all
    final_output = [rule_dict[r] for r in opt_rules if config['rules'][r]]

    if NCBI == []:
        pass
    else:
        final_output.extend(expand("data/processed/{name}/tables/df_ncbi_meta.csv", name = PROJECT_IDS))
    return final_output

def get_final_output(df_samples):
    """
    Generate outputs of for all projects
    """
    sys.stderr.write(f"Step 6. Preparing list of final outputs...\n")
    final_output = []
    for p in config["projects"]:
        sys.stderr.write(f" - Getting outputs for project: {p['name']}\n")
        try:
            final_output.extend(get_project_outputs(p['rules'], p['name'], df_samples))
            sys.stderr.write(f"   {p['name']}: Using project rule configuration.\n")
        except KeyError:
            final_output.extend(get_project_outputs(config, p['name'], df_samples))
            sys.stderr.write(f"   {p['name']}: Using global rule configuration.\n")
    sys.stderr.write(f" - Ready to generate all outputs.\n\n")
    return final_output

##### 8. Set up custom resource directory provided in config["resources_path"] #####
def custom_resource_dir():
    """
    Generate symlink for user defined resources location
    """
    resource_dbs = config["resources_path"]
    sys.stderr.write(f"Step 5. Checking for user-defined local resources...\n")
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
                    sys.stderr.write(f" - Symlink for {r} already generated from: {existing_path}\n")
                # update symlink because new path is given
                else:
                    slink.unlink()
                    slink.symlink_to( path )
                    updated_path = Path.readlink(Path(slink)) 
                    sys.stderr.write(f" - Updating symlink for {r} from: {existing_path} to: {updated_path}\n")
            # generate a new symlink
            except FileNotFoundError:
                sys.stderr.write(f" - Generating symlink for {r} from: {path}\n")                    
                slink.symlink_to( path )
        # raise an Error if external path not found
        else:
            raise FileNotFoundError(f"Error: User-defined resource {r} at {path} does not exist.\nCheck the config.yaml and provide the right path for resource {r} or\nchange it to the default path: resources/{r}\n")
    sys.stderr.write(f"   All resources set.\n\n")
    return
