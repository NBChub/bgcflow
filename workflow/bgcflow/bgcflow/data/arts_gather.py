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


def combine_arts_json(input_json, change_log_path, table):
    logging.info("Combining and correcting ARTS output...")

    input_json = Path(input_json)
    logging.info(input_json)
    if input_json.is_file() and input_json.suffix == ".json":
        logging.info(f"Getting ARTS overview from a single file: {input_json}")
        input_json_files = input_json

    elif input_json.is_file() and input_json.suffix == ".txt":
        logging.info(f"Getting ARTS overview  from a text file: {input_json}")
        with open(input_json, "r") as file:
            file_content = [i.strip("\n") for i in file.readlines()]
            if len(file_content) == 1:
                # Paths space-separated on a single line
                paths = file_content[0].split()
            else:
                # Paths written on separate lines
                paths = file_content
            input_json_files = [
                Path(path) for path in paths if Path(path).suffix == ".json"
            ]
    else:
        input_json_files = [
            Path(file)
            for file in str(input_json).split()
            if Path(file).suffix == ".json"
        ]
        logging.info(
            f"Getting ARTS overview from the given list of {len(input_json_files)} files..."
        )

    container = {}
    change_dict = generate_change_dict(change_log_path)

    with alive_bar(len(input_json_files), title="Merging json:") as bar:
        for j in input_json_files:
            arts_json = Path(j)
            with open(arts_json, "r") as f:
                value = json.load(f)
            container.update(value)
            bar()

    df = pd.DataFrame.from_dict(container).T
    df.index.name = "pkey"

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(table)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    combine_arts_json(sys.argv[1], sys.argv[2], sys.argv[3])
