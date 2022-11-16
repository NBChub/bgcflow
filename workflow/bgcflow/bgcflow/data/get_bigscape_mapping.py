from pathlib import Path
import pandas as pd
import sys


def get_bigscape_mapping(path, outfile):
    """
    Given a path containing directories of antiSMASH result, return a table
    which maps BGC regions with its parent folder (aka accession ids)
    """
    container = {}
    for x in Path(path).iterdir():
        if x.is_dir():
            genome_id = x.name
            bgcs = {item.stem: genome_id for item in list(x.glob("*.region*.gbk"))}
            container.update(bgcs)
    df = pd.DataFrame.from_dict(container, orient="index", columns=["genome_id"])
    df.index.name = "bgc_id"
    df.to_csv(outfile)
    return None


if __name__ == "__main__":
    get_bigscape_mapping(sys.argv[1], sys.argv[2])
