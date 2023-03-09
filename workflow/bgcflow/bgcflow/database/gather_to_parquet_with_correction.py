import json
import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def correct_bgc_id_overview(
    overview_file, mapping_file, genome_id=False, exclude="_bgc_overview.json"
):
    """
    Use a mapping file to correct bgc_ids
    """
    logging.info(f"Correcting shortened bgc ids for {genome_id}...")
    overview_path = Path(overview_file)
    mapping_path = Path(mapping_file)

    with open(overview_path, "r") as f:
        overview = json.load(f)

    with open(mapping_path, "r") as f:
        mapping = json.load(f)

    if not genome_id:
        genome_id = overview_path.strip(exclude)
    else:
        pass

    new_dict = {}

    # correct index
    for keys in overview.keys():
        for m in mapping[genome_id].keys():
            if keys in m:
                log_change = mapping[genome_id][m]
                correct_bgc_id = Path(log_change["symlink_path"]).stem
                if log_change["record_id"] != log_change["original_id"]:
                    logging.debug(f"Replacing Index: {keys} to {correct_bgc_id}")
                    new_dict[correct_bgc_id] = overview[keys]
                    new_dict[correct_bgc_id]["accession"] = log_change["original_id"]
                else:
                    new_dict[correct_bgc_id] = overview[keys]
                    pass

        for col in overview[keys].keys():
            if col in ["region_id", "bgc_id"]:
                value = overview[keys][col]
                if value is not None:
                    for m in mapping[genome_id].keys():
                        if value in m:
                            log_change = mapping[genome_id][m]
                            correct_bgc_id = Path(log_change["symlink_path"]).stem
                            if log_change["record_id"] != log_change["original_id"]:
                                logging.debug(
                                    f"Replacing {keys}[{col}]: {value} to {correct_bgc_id}"
                                )
                                overview[keys][col] = correct_bgc_id
                                new_dict[keys] = overview[keys]
                            else:
                                new_dict[keys] = overview[keys]
                                pass

    logging.info("Done!")
    return new_dict


def gather_to_parquet(
    input_json, mapping_dir, table, exclude="_bgc_overview.json", index_name="bgc_id"
):
    input_json = input_json.split()

    merged_dict = {}
    for j in input_json:
        mapping_file = Path(j)
        genome_id = mapping_file.name.replace(exclude, "")
        mapping_path = Path(mapping_dir) / f"{genome_id}/{genome_id}-change_log.json"
        corrected = correct_bgc_id_overview(mapping_file, mapping_path, genome_id)
        merged_dict.update(corrected)

    df = pd.DataFrame.from_dict(merged_dict).T
    df.index.name = index_name

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(table)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    gather_to_parquet(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
