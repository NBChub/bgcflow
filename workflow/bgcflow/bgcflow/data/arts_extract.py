import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def convert_bgc_id(name, genome_id):
    cluster_number = name.strip("cluster-").split("_")
    bgc_id = f"{genome_id}.region{str(cluster_number[-1]).zfill(3)}"
    return bgc_id


def arts_bgc_hits(tsv, genome_id, outfile):
    logging.debug(f"Extracting ARTS2 result for {genome_id}")
    tsv = Path(tsv)
    df = pd.read_csv(tsv, sep="\t")
    df["#Cluster"] = [convert_bgc_id(c, genome_id) for c in df["#Cluster"]]
    df.rename(columns={"#Cluster": "bgc_id"}).set_index("bgc_id").T.to_json(outfile)
    return


if __name__ == "__main__":
    arts_bgc_hits(sys.argv[1], sys.argv[2], sys.argv[3])
