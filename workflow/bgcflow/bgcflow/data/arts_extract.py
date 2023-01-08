import json
import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def generate_arts_dict(df_arts):
    """
    Convert arts table to dictionary
    """
    arts_dict = {}
    arts_scaffold_count = df_arts.Source.value_counts().to_dict()
    for k in arts_scaffold_count.keys():
        arts_dict[k] = {"counts": arts_scaffold_count[k]}
        arts_dict[k]["regions"] = []
        for i in df_arts[df_arts.loc[:, "Source"] == k].index:
            regions = {
                "cluster_id": df_arts.loc[i, "#Cluster"],
                "products": df_arts.loc[i, "Type"],
                "location": df_arts.loc[i, "Location"],
            }
            arts_dict[k]["regions"].append(regions)
    return arts_dict


def generate_arts_mapping(df_arts, as_json):
    arts_dict = generate_arts_dict(df_arts)
    # container for final result
    hit_mapper = {}
    contig_mapper = {}

    # iterate antismash json records
    for num, r in enumerate(as_json["records"]):
        # count number of detected regions per record
        region_count = len(r["areas"])
        # if antismash detects region
        if region_count > 0:
            contig_id = r["id"]
            logging.info(
                f"Contig {contig_id} has {region_count} regions detected. Finding corresponding scaffold in ARTS2..."
            )
            # find arts scaffold with the same region number detected
            arts_match = [
                k
                for k in arts_dict.keys()
                if int(arts_dict[k]["counts"]) == int(region_count)
            ]
            logging.debug(f"Finding matches from: {arts_match}")
            for n, a in enumerate(r["areas"]):
                bgc_id = f"{contig_id}.region{str(n+1).zfill(3)}"
                location = f"{a['start']} - {a['end']}"
                products = ",".join(a["products"])
                logging.debug(f"Finding match for: {bgc_id} | {location} | {products}")
                for m in arts_match:
                    bgc_match = arts_dict[m]["regions"][n]
                    if (bgc_match["location"] == location) and (
                        bgc_match["products"] == products
                    ):
                        logging.debug(
                            f"Found match! {bgc_match['cluster_id']} | {bgc_match['location']} | {bgc_match['products']}"
                        )
                        # logging.debug(f"Found match! {contig_id} == ARTS2 {m}")
                        hit_mapper[bgc_match["cluster_id"]] = bgc_id
                        contig_mapper[bgc_match["cluster_id"]] = contig_id
    return hit_mapper, arts_dict, contig_mapper


def extract_arts_table(input_table, input_as_json, outfile, genome_id=None):
    """
    Given an ARTS2 bgc table, map the corresponding cluster hits to antismash regions
    """
    df_arts = pd.read_csv(input_table, sep="\t")

    with open(input_as_json, "r") as f:
        as_json = json.load(f)

    if genome_id is None:
        logging.info("Assuming genome_id from ARTS2 input")
        genome_id = Path(input_table).parents[1].stem

    logging.info(f"Extracting ARTS2 BGC hits from: {genome_id}")

    logging.info("Generating ARTS2 to antiSMASH region mapper...")
    hit_mapper, arts_dict, contig_mapper = generate_arts_mapping(df_arts, as_json)

    logging.info("Mapping ARTS2 cluster to antiSMASH regions...")
    df = pd.DataFrame(columns=df_arts.columns)
    for i in df_arts.index:
        # only extract hits
        if len(df_arts.loc[i, "Genelist"]) > 2:
            cluster_id = df_arts.loc[i, "#Cluster"]
            bgc_id = hit_mapper[cluster_id]
            contig_id = contig_mapper[cluster_id]
            df.loc[i, :] = df_arts.loc[i, :]
            df.loc[i, "#Cluster"] = bgc_id
            df.loc[i, "Source"] = contig_id
            df.loc[i, "genome_id"] = genome_id
    df = df.rename(columns={"#Cluster": "bgc_id"}).set_index("bgc_id")

    logging.info(f"Writing results to {outfile}")
    df.T.to_json(outfile)
    return


if __name__ == "__main__":
    extract_arts_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
