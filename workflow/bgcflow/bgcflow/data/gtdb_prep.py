import os, sys, json
import pandas as pd
import requests

def gtdb_prep(genome_id, outfile, df_path, release='R202'):
    """
    Given a genome id and the samples Pandas dataframe, write  a JSON file containing taxonomic information from GTDB API. The script will first search using the closest taxonomic placement (NCBI accession id), then using the genus information provided by user. If no information is provided, return an empty taxonomic information.
    """
    df = pd.read_csv(df_path)
    query = df[df.loc[:, "genome_id"] == genome_id].fillna("")

    gtdb_tax = {}
    sys.stderr.write(f"Fetching taxonomic information for {genome_id}...\n")
    if query.source.values[0] == "ncbi":
        gtdb_tax = get_ncbi_taxon_GTDB(query.genome_id.values[0], release)
    elif query.closest_placement_reference.values[0] != "":
        gtdb_tax = get_ncbi_taxon_GTDB(query.closest_placement_reference.values[0], release)
        gtdb_tax["genome_id"] = genome_id
    elif query.genus.values[0] != "":
        gtdb_tax['genome_id'] = genome_id
        gtdb_tax.update(get_parent_taxon_GTDB(query.genus.values[0], "genus", release))
        gtdb_tax['gtdb_taxonomy']["species"] = f"s__{gtdb_tax['gtdb_taxonomy']['genus'].split('__')[-1]} sp."
    else:
        sys.stderr.write(f"No taxonomic information found, returning an empty dictionary for {genome_id}.\n")
        gtdb_tax['genome_id'] = genome_id
        gtdb_tax["gtdb_taxonomy"] = {"domain" : "",
                                     "phylum" : "",
                                     "class" : "",
                                     "order" : "",
                                     "family" : "",
                                     "genus" : "",
                                     "species" : ""}

    with open(outfile, 'w') as file:
        json.dump(gtdb_tax, file, indent=2)
        file.close
    return

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
    gtdb_prep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])