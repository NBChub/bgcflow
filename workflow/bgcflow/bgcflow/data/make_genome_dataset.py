import json
import logging
import sys
from pathlib import Path

import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def combine_bgc_counts(input_json, filter_str="_bgc_counts"):
    container = {}
    logging.info("Reading json files...")
    with alive_bar(len(input_json), title="Merging BGC counts:") as bar:
        for item in input_json:
            item = Path(item)
            genome_id = item.stem
            if filter_str in genome_id:
                genome_id = genome_id.replace(filter_str, "")
            logging.debug(f"Processing {genome_id}")
            with open(item, "r") as f:
                data = json.load(f)
                products = data[genome_id].pop("products")
                data[genome_id] = data[genome_id] | products
                container.update(data)
            bar()
    return container


def write_genome_table(input_json, samples_table, genome_table):
    """
    Write df_genomes.csv table in processed data
    """
    # Accomodate multiple inputs to generate dataframe
    shell_input = samples_table.split()
    logging.info(f"Reading samples table: {shell_input}")
    dfList = [pd.read_csv(s).set_index("genome_id", drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)

    # Handle multiple json
    input_json = Path(input_json)
    logging.info(input_json)
    if input_json.is_file() and input_json.suffix == ".json":
        logging.info(f"Getting BGC overview from a single file: {input_json}")
        input_json_files = input_json

    elif input_json.is_file() and input_json.suffix == ".txt":
        logging.info(f"Getting BGC overview  from a text file: {input_json}")
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
            f"Getting BGC overview from the given list of {len(input_json_files)} files..."
        )

    bgc_counts = combine_bgc_counts(input_json_files)
    bgc_counts = pd.DataFrame.from_dict(bgc_counts).T

    logging.debug(f"Writing file to: {genome_table}")

    # Generate dataframe
    df_genomes = pd.concat([df_samples, bgc_counts], axis=1)

    # Save dataframes to csv tables
    genome_table = Path(genome_table)
    genome_table.parent.mkdir(parents=True, exist_ok=True)
    df_genomes.to_csv(genome_table, index=False)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    write_genome_table(sys.argv[1], sys.argv[2], sys.argv[3])
