import logging
import sys
from pathlib import Path

import networkx as nx
import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def extract_cog_information(mmseqs2_cluster_tsv, outfile):
    """
    Extract COG information from mmseqs2
    """
    G = nx.read_edgelist(Path(mmseqs2_cluster_tsv))
    # get the connected components of G
    components = nx.connected_components(G)

    # sort the components based on their size
    sorted_components = sorted(components, key=lambda x: len(x), reverse=True)

    result = {"feature_id": [], "cluster_id": [], "cluster_n": [], "locus_tag": []}
    for num, g in enumerate(sorted_components):
        size = len(g)
        for item in g:
            result["feature_id"].append(f"cds-{item}")
            result["cluster_id"].append(f"cog_{num+1:02d}")
            result["cluster_n"].append(size)
            result["locus_tag"].append(item)
    df_mmseqs2 = pd.DataFrame.from_dict(result)
    df_mmseqs2.to_csv(outfile, index=False)
    return


if __name__ == "__main__":
    extract_cog_information(sys.argv[1], sys.argv[2])
