import pandas as pd
from pathlib import Path
import sqlite3
import logging
from alive_progress import alive_bar
import json, sys

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def gcf_hits(df_gcf, cutoff=900):
    """
    Filter bigslice result based on distance threshold to model and generate data.
    """
    mask = df_gcf.loc[:, "membership_value"] <= cutoff
    df_gcf_filtered = df_gcf[mask]
    bgcs = df_gcf_filtered.bgc_id.unique()
    gcfs = df_gcf_filtered.gcf_id.unique()
    logging.debug(f"BiG-SLICE query with BiG-FAM run 6, distance cutoff {cutoff}")
    logging.debug(f"Number of bgc hits : {len(bgcs)}/{len(df_gcf.bgc_id.unique())}")
    logging.debug(f"Number of GCF hits : {len(gcfs)}")
    return df_gcf_filtered

def grab_gcf_model_summary(gcf_query, database):
    """
    Summarize gcf model from a give list of bigfam gcf ids
    gcf_query = list of bigfam gfc id
    database = the full run database of bigslice
    """
    
    def grab_mibig_members(q, mibig_bigfam):
        """
        Summarize information of mibig members in a model
        """
        mibig_members = list(set(q.bgc_id) & set(mibig_bigfam.id))
        mibig_members = [mibig_bigfam.loc[i, "name"] for i in mibig_members]
        return mibig_members
    
    # initiate connection
    logging.info(f"Initiating SQL connection to {database}...")
    conn = sqlite3.connect(database)

    # load main tables
    logging.info(f"Reading gcf_membership table...")
    df_bigfam = pd.read_sql_query(f"select * from gcf_membership;", conn)
    logging.info(f"Reading bgc table...")
    df_bgc_bigfam = pd.read_sql_query(f"select * from bgc;", conn)

    # Load BiG-FAM dataset
    mibig_bigfam = df_bgc_bigfam[df_bgc_bigfam.loc[:, "type"] == "mibig"].set_index("id", drop=False)

    # Filter for hit queries
    df_bigfam_query = df_bigfam[df_bigfam.loc[:, 'gcf_id'].isin(gcf_query)]

    # return information of each model
    logging.info(f"Summarizing information for {len(gcf_query)} GCF models...")
    summary = {}
    with alive_bar(len(gcf_query)) as bar:
        for g in gcf_query:
            values = {}

            q = df_bigfam[df_bigfam.loc[:, 'gcf_id'] == g]
            q = q[q.loc[:, 'rank'] == 0] # select only first hits

            # get core member info
            q_core = q[q.loc[:, 'membership_value'] <= 900]
            q_core_mibig = grab_mibig_members(q_core, mibig_bigfam)

            # get putative member info
            q_putative = q[q.loc[:, 'membership_value'] > 900]
            q_putative_mibig = grab_mibig_members(q_putative, mibig_bigfam)

            values["core_member"] = len(q_core)
            values["putative_member"] = len(q_putative)
            values["core_member_mibig"] = q_core_mibig
            values["putative_member_mibig"] = q_putative_mibig
            values["core_member_mibig_count"] = len(q_core_mibig)
            values["putative_member_mibig_count"] = len(q_putative_mibig)
            values["link to BiG-FAM"] = f"https://bigfam.bioinformatics.nl/run/6/gcf/{g}"
            summary[str(g)] = values
            
            bar()
    
    return summary

def summarize_bigslice_query(bigslice_query_path, output_path, database_path="resources/bigslice/full_run_result/result/data.db", cutoff=900):
    """
    Summarize bigslice query result.
    
    Input:
    - bigslice_query_path
    
    Output:
    - A folder to store two output files:
        1. query_network.csv to build Cytoscape Network
        2. JSON file summarizing GCF model hits
    """
    bigslice_query = Path(bigslice_query_path)
    database = Path(database_path)
    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)
    
    df_gcf_membership = pd.read_csv(bigslice_query / "gcf_membership.csv")
    df_bgc = pd.read_csv(bigslice_query / "bgc.csv")

    # Create edge table network for cytoscape enriched with ids for mapping metadata (bgc_name, genome_id)
    logging.info("Generating network table of BiG-SLICE hits...")
    data = gcf_hits(df_gcf_membership, cutoff)
    bgc_info = df_bgc.set_index("id", drop=True)
    bgc_info["bgc_id"] = [str(i).replace(".gbk", "") for i in bgc_info["orig_filename"]]
    for i in data.index:
        bgc_id = data.loc[i, "bgc_id"]
        data.loc[i, "bgc_id"] = str(bgc_info.loc[bgc_id, "bgc_id"])
    data.to_csv(output / "query_network.csv")
    
    # Summarizing GCF model hits
    logging.info("Summarizing GCF model hits...")
    gcf_query = list(data.gcf_id.unique())
    gcf_summary = grab_gcf_model_summary(gcf_query, database_path)

    with open(output / "gcf_summary.json", "w") as outfile:
        json.dump(gcf_summary, outfile, indent=4)
    
    logging.info("Job done")
    return

if __name__ == "__main__":
    summarize_bigslice_query(sys.argv[1], sys.argv[2])