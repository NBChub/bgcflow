import json
import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def correct_bgc_id_overview(overview_file, mapping_file, genome_id=False):
    """
    Use a mapping file to correct bgc_ids

    Parameters
    ----------
    overview_file : str / path
        Path to the BGC overview file in JSON format

    mapping_file : str / path
        Path to the mapping file in JSON format

    genome_id : str, optional
        Genome ID to correct BGC IDs for, if not provided, it is extracted from the overview file name

    Returns
    -------
    new_dict : dict
        Corrected BGC overview dictionary with updated BGC IDs
    """
    logging.info(f"Correcting shortened bgc ids for {genome_id}...")
    overview_path = Path(overview_file)
    mapping_path = Path(mapping_file)

    with open(overview_path, "r") as f:
        overview = json.load(f)

    with open(mapping_path, "r") as f:
        mapping = json.load(f)

    if not genome_id:
        genome_id = overview_path.stem.strip("_bgc_overview")
    else:
        pass

    new_dict = {}

    for bgc_id in overview.keys():
        for m in mapping[genome_id].keys():
            if bgc_id in m:
                log_change = mapping[genome_id][m]
                correct_bgc_id = Path(log_change["symlink_path"]).stem
                if log_change["record_id"] != log_change["original_id"]:
                    logging.debug(f"Replacing {bgc_id} to {correct_bgc_id}")
                    new_dict[correct_bgc_id] = overview[bgc_id]
                    new_dict[correct_bgc_id]["accession"] = log_change["original_id"]
                else:
                    new_dict[correct_bgc_id] = overview[bgc_id]
                    pass
                new_dict[correct_bgc_id]["source"] = "bgcflow"
                new_dict[correct_bgc_id]["gbk_path"] = log_change["target_path"]
    logging.info("Done!")
    return new_dict


def gather_bgc_overview(input_json, mapping_dir, table):
    """
    Gather BGC overview data from multiple JSON files and create a merged table

    Parameters
    ----------
    input_json : str
        Two different input types can be used:
        - Space-separated paths to BGC overview JSON files (put the string inside '' in bash expression)
        - A text file (.txt) containing the paths of the json files (one per line, or space-separated on a single line)
    mapping_dir : str / path
        Directory containing mapping files

    table : str / path
        Path to the output merged table

    Returns
    -------
    None
    """
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

    merged_dict = {}
    for j in input_json_files:
        mapping_file = Path(j)
        genome_id = mapping_file.name.replace("_bgc_overview.json", "")
        mapping_path = Path(mapping_dir) / f"{genome_id}/{genome_id}-change_log.json"
        corrected = correct_bgc_id_overview(mapping_file, mapping_path, genome_id)
        merged_dict.update(corrected)

    df = pd.DataFrame.from_dict(merged_dict).T
    df.index.name = "bgc_id"

    logging.debug("Checking similarity values...")
    df["similarity"] = df["similarity"].apply(lambda x: 1 if x > 1 else x)

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(table)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    gather_bgc_overview(sys.argv[1], sys.argv[2], sys.argv[3])
