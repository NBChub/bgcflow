import json
import logging
import sys
from pathlib import Path

import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def generate_change_dict(change_log):
    logging.info("Generate BGC ids mapping dict...")
    change_log = Path(change_log)
    change_log = [i for i in change_log.glob("*/*.json")]
    change_dict = {}
    for i in change_log:
        with open(i, "r") as f:
            data = json.load(f)
            for k in data.keys():
                change_dict[k] = data[k]
    return change_dict


def correct_arts_bgc_ids(arts_json, change_dict, genome_id):
    logging.debug(f"Correcting BGC ids for {genome_id}")
    output = {}
    with open(arts_json, "r") as f:
        query = json.load(f)
        for k in query.keys():
            correct_bgc_id = Path(
                change_dict[genome_id][f"{k}.gbk"]["symlink_path"]
            ).stem
            output[correct_bgc_id] = query[k]
    return output


def combine_arts_json(input_json, change_log_path, table):
    logging.info("Combining and correcting ARTS output...")

    # Handle multiple json
    input_json = input_json.split()

    container = {}
    change_dict = generate_change_dict(change_log_path)

    with alive_bar(len(input_json), title="Merging json:") as bar:
        for j in input_json:
            arts_json = Path(j)
            genome_id = arts_json.stem
            value = correct_arts_bgc_ids(arts_json, change_dict, genome_id)
            container.update(value)
            bar()

    df = pd.DataFrame.from_dict(container).T
    df.index.name = "bgc_id"

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(table)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    combine_arts_json(sys.argv[1], sys.argv[2], sys.argv[3])
