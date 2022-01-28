import os, sys, json
import numpy as nan
import pandas as pd
import requests

def gtdb_prep(genome_id, outfile, samples_table, tax_path, release='R202'): # what happen if it does not find?
    """
    Given a genome id and the samples Pandas dataframe, write  a JSON file containing taxonomic information from GTDB API. The script will first search using the closest taxonomic placement (NCBI accession id), then using the genus information provided by user. If no information is provided, return an empty taxonomic information.
    """

    class EmptyTaxError(Exception):
        """Raised when this script returns empty dict"""
        pass

    class PlacementError(Exception):
        """Raised when this script returns empty dict"""
        pass

    def empty_result(genome_id):
        """
        Helper script for creating empty result
        """
        sys.stderr.write(f"No taxonomic information found, returning empty values for {genome_id}.\n")
        gtdb_tax = {}
        gtdb_tax['genome_id'] = genome_id
        gtdb_tax["gtdb_taxonomy"] = {"domain" : "d__",
                                    "phylum" : "p__",
                                    "class" : "c__",
                                    "order" : "o__",
                                    "family" : "f__",
                                    "genus" : "g__",
                                    "species" : "s__"}
        return gtdb_tax


    def find_taxonomy(query, genome_id, gtdb_tax):
        """
        Helper script to decide taxonomic placement for a given query
        """
        # If closest placement reference is provided, try finding taxonomic information from GTDB API
        if query.closest_placement_reference.values[0] != "":
            try:
                sys.stderr.write("Inferring taxonomic placement from provided closest reference....\n")
                gtdb_tax = get_ncbi_taxon_GTDB(query.closest_placement_reference.values[0], release)
                gtdb_tax["genome_id"] = genome_id
            except KeyError as e:
                raise PlacementError(f"Cannot infer taxonomic placement from provided closest reference. Make sure the accession id: {query.closest_placement_reference.values[0]} is part of GTDB release: {e}\n")


        # If NCBI accession is provided, try to find taxonomic information from GTDB API
        elif query.source.values[0] == "ncbi":
            try:
                sys.stderr.write("Inferring taxonomic placement from NCBI accession....\n")
                gtdb_tax = get_ncbi_taxon_GTDB(query.genome_id.values[0], release)
            except KeyError:
                if query.genus.values[0] != "":
                    sys.stderr.write("Inferring taxonomic placement from provided genus information....\n")
                    gtdb_tax['genome_id'] = genome_id
                    gtdb_tax.update(get_parent_taxon_GTDB(query.genus.values[0], "genus", release))
                    gtdb_tax['gtdb_taxonomy']["species"] = f"s__{gtdb_tax['gtdb_taxonomy']['genus'].split('__')[-1]} sp."
                else:
                    gtdb_tax = empty_result(genome_id)

        # Try to get taxonomic information from genus information
        elif query.genus.values[0] != "":
            sys.stderr.write("Inferring taxonomic placement from provided genus information....\n")
            gtdb_tax['genome_id'] = genome_id
            gtdb_tax.update(get_parent_taxon_GTDB(query.genus.values[0], "genus", release))
            gtdb_tax['gtdb_taxonomy']["species"] = f"s__{gtdb_tax['gtdb_taxonomy']['genus'].split('__')[-1]} sp."

        # If no information is found, return an empty dict
        else:
            gtdb_tax = empty_result(genome_id)

        return gtdb_tax

    # get query by subsetting samples df with genome id
    shell_input = samples_table.split()
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)
    query = df_samples[df_samples.loc[:, "genome_id"] == genome_id].fillna("")

    # create empty container
    gtdb_tax = {}

    # Starting process
    sys.stderr.write(f"Fetching taxonomic information for {genome_id}...\n")
    # Go through user provided taxonomic placement
    if any(os.path.isfile(t) for t in tax_path.split()):
        try:
            gtdb_tax = get_user_defined_classification(genome_id, tax_path)
            print(gtdb_tax)
        except KeyError:
            sys.stderr.write(f"{genome_id}: Not found in user provided taxonomic placement...\n")
            gtdb_tax = find_taxonomy(query, genome_id, gtdb_tax)
            print(gtdb_tax)
    else:
        gtdb_tax = find_taxonomy(query, genome_id, gtdb_tax)

    if gtdb_tax == {}:
        raise EmptyTaxError(f"Oops, this shouldn't happen. It returns an empty dict. Something is wrong with the script.")

    with open(outfile, 'w') as file:
        json.dump(gtdb_tax, file, indent=2)
        file.close
    return

def get_user_defined_classification(genome_id, tax_path):
    """
    Get taxonomic information from user provided GTDB-like output
    """
    shell_input = tax_path.split()
    dfList = [pd.read_csv(s, sep="\t").set_index('user_genome', drop=False) for s in shell_input]
    df_tax = pd.concat(dfList, axis=0)

    level_dict = {"d" :"domain",
                  "p" : "phylum",
                  "c" : "class",
                  "o" : "order",
                  "f" : "family",
                  "g" : "genus",
                  "s" : "species"}

    query = df_tax.loc[genome_id, "classification"].split(";")

    result = {}
    result["genome_id"] = genome_id
    result["gtdb_url"] = "user provided classification"
    result["gtdb_release"] = "unknown"
    result["gtdb_taxonomy"] = {level_dict[q.split("__")[0]] : q for q in query}
    sys.stderr.write(f"Using user provided GTDB classification.\n")
    return result

def get_ncbi_taxon_GTDB(accession, release = "R202"):
    """
    Given an NCBI accession, return a json object of taxonomic information from GTDB API
    """
    api_url = f"https://gtdb.ecogenomic.org/api/v1/taxonomy/genome/{accession}"
    response = requests.get(api_url)
    js = response.json()
    result = {}
    result["genome_id"] = accession
    result["gtdb_url"] = api_url
    result["gtdb_release"] = release
    try:
        result["gtdb_taxonomy"] = js["gtdb_taxonomy"][release]
    except KeyError as err:
        if err.args[0] == "gtdb_taxonomy":
            sys.stderr.write(f"Malformed genome id: {accession}. Make sure to use the right NCBI genome accession format.\n")
            raise
        elif err.args[0] == release:
            sys.stderr.write(f"Cannot find genome id: {accession} in GTDB API.\n")
            raise
    return result

def get_parent_taxon_GTDB(taxon, level, release = "R202"):
    """
    Given a taxon and its level, return a json object of parent taxons from GTDB API
    """
    level_dict = {"domain" : "d",
                  "phylum" : "p",
                  "class" : "c",
                  "order" : "o",
                  "family" : "f",
                  "genus" : "g",
                  "species" : "s"}

    try:
        query = f"{level_dict[level]}__{taxon}"
    except KeyError:
        sys.stderr.write(f"Incorrect taxon level format. Please choose from available format: {list(level_dict.keys())}.\n")
        raise

    api_url = f"https://gtdb.ecogenomic.org/api/v1/taxonomy/taxon/{query}/parent"
    response = requests.get(api_url)
    js = response.json()
    result = {}
    result["gtdb_url"] = api_url
    result["gtdb_release"] = release
    result["gtdb_taxonomy"] = {}

    try:
        for k in js.keys():
            if release in js[k]["releases"]:
                t = js[k]["taxon"]
                result["gtdb_taxonomy"][k] = t
    except TypeError:
        sys.stderr.write(f"{js['message']}\n")
        raise

    return result

if __name__ == "__main__":
    gtdb_prep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])