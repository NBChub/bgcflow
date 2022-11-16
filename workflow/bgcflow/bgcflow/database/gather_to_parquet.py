from pathlib import Path
import sys
import pandas as pd
from alive_progress import alive_bar
import json
import logging

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def combine_json(input_json):
    container = {}
    logging.info(f"Reading json files...")
    with alive_bar(len(input_json), title="Merging json:") as bar:
        for item in input_json:
            item = Path(item)
            logging.debug(f"Processing {item.stem}")
            with open(item, "r") as f:
                data = json.load(f)
                container.update(data)
            bar()
    return container


def write_parquet(input_json, index_key, table):
    """
    Write .parquet table in processed data
    """
    # Handle multiple json
    input_json = input_json.split()
    df = combine_json(input_json)
    df = pd.DataFrame.from_dict(df).T
    df.index.name = index_key

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(table)
    logging.info(f"Job done")
    return None


if __name__ == "__main__":
    write_parquet(sys.argv[1], sys.argv[2], sys.argv[3])
