import logging
import sys
from pathlib import Path

import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def generate_bgc_id_mapping(changelog, genome_id, contig=1):
    df_changelog = pd.read_csv(changelog)
    df_mapping = df_changelog[df_changelog.genome_id == genome_id]
    mapping = {}
    region_ctr = 1
    for bgc_id in df_mapping.bgc_id:
        region = int(bgc_id.split("region")[-1])
        if region == region_ctr:
            mapping[f"cluster-{contig}_{region}"] = bgc_id
            region_ctr = region_ctr + 1
            pass
        else:
            region_ctr = 1  # reset
            contig = contig + 1
            mapping[f"cluster-{contig}_{region}"] = bgc_id
            region_ctr = region_ctr + 1
    return mapping


def extract_arts2(input_folder, changelog, outfile):
    # Output handler for aggregation
    df_hits_container = []

    # Handle multiple folders
    input_arts = {Path(i) for i in input_folder.split()}

    with alive_bar(len(input_arts), title="Extracting ARTS2:") as bar:
        for arts in input_arts:
            genome_id = arts.name
            logging.info(f"Extracting ARTS2 result from {genome_id}")
            tsv = arts / "tables/bgctable.tsv"
            df_hits = pd.read_csv(tsv, sep="\t")
            mapping = generate_bgc_id_mapping(changelog, genome_id)
            df_hits = df_hits.rename(columns={"#Cluster": "bgc_id"})
            df_hits["genome_id"] = genome_id
            try:
                for i in df_hits.index:
                    c = df_hits.loc[i, "bgc_id"]
                    logging.debug(f"{c}, {mapping[c]}")
                    df_hits.loc[i, "bgc_id"] = mapping[c]
            except KeyError:
                logging.warning(f"{c} didn't found any matches")
                logging.debug(f"Available keys: {list(mapping.keys())}")
                assert len(df_hits.bgc_id) == len(mapping)
                logging.debug(
                    "ARTS2 ids and mapping has the same length. Adjusting mapping keys..."
                )
                starting_contig = int(c.strip("cluster-").split("_")[0])
                new_mapping = generate_bgc_id_mapping(
                    changelog, genome_id, contig=starting_contig
                )
                logging.debug("Rerunning extraction...")
                for i in df_hits.index:
                    c = df_hits.loc[i, "bgc_id"]
                    logging.debug(f"{c}, {new_mapping[c]}")
                    df_hits.loc[i, "bgc_id"] = new_mapping[c]
            df_hits_container.append(df_hits)
            bar()

    logging.info(f"Merging all data into final output: {outfile}")
    pd.concat(df_hits_container).set_index("bgc_id").to_csv(outfile)
    return


if __name__ == "__main__":
    extract_arts2(sys.argv[1], sys.argv[2], sys.argv[3])
