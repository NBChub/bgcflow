import os
import sys
import pandas as pd
import requests

def get_parent_taxon_GTDB(taxon, level, release = "R202"):
    """
    Given a taxon and its level, return a json object of parent taxons from GTDB API
    """
    level_dict = {"genus" : "g",
                  "species" : "s"}
    query = f"{level_dict[level]}__{taxon}"
    api_url = f"https://gtdb.ecogenomic.org/api/v1/taxonomy/taxon/{query}/parent"
    response = requests.get(api_url)
    js = response.json()
    result = {}
    for k in js.keys():
        if release in js[k]["releases"]:
            t = js[k]["taxon"]
            result[k.title()] = t.split("__")[-1]
    return result

def fix_gtdb_taxonomy(samples_path, taxonomy_raw, taxonomy, meta):
    shell_input = samples_path.replace("[","").replace("]","")
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input.split()]
    df_samples = pd.concat(dfList, axis=0).fillna(0)

    df_tax = pd.read_csv(taxonomy_raw, sep="\t").fillna("")
    df_tax = df_tax.set_index("#Genome folder", drop=False)

    dict_tax = {i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("custom")].genome_id}
    dict_tax.update({i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("patric")].genome_id})
    dict_tax.update({i : str(i)+"/" for i in df_samples[df_samples.source.eq("ncbi")].genome_id})
    container = []
    df_meta = pd.DataFrame(index=df_samples.index, columns=df_tax.columns)
    for i in df_samples.genome_id:
        try:
            line = df_tax.loc[dict_tax[i], :].to_dict()
            df_meta.loc[i, :] = df_tax.loc[dict_tax[i], :]
        except KeyError:
            line = {'#Genome folder' : "",
                    'Kingdom' : "",
                    'Phylum' : "",
                    'Class' : "",
                    'Order' : "",
                    'Family' : "",
                    'Genus' : "",
                    'Species' : "",
                    'Organism' : ""
                    }
            df_meta.loc[i, :] = line

        line.update({'#Genome folder': f"{i}/"})

        # Try to fix with value from genus
        if df_meta.loc[i, "Genus"] == "" and df_samples.loc[i, "genus"] != 0:
            line = get_parent_taxon_GTDB(df_samples.loc[i, "genus"], "genus")
            line["Kingdom"] = line.pop("Domain")
            line["Species"] = "sp."
            line["Organism"] = f"{line['Genus']} {line['Species']}"
            line["#Genome folder"] = f"{i}/"
            df_meta.loc[i, :] = line

        container.append(line)

    df = pd.DataFrame(container)
    df.to_csv(taxonomy, sep="\t", index=False)
    df_meta.drop(columns='#Genome folder', inplace=True)
    df_meta.to_csv(meta)
    return None

if __name__ == "__main__":
    fix_gtdb_taxonomy(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])