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
    input_json = input_json.split()

    merged_dict = {}
    for j in input_json:
        mapping_file = Path(j)
        genome_id = mapping_file.name.replace("_bgc_overview.json", "")
        mapping_path = Path(mapping_dir) / f"{genome_id}/{genome_id}-change_log.json"
        corrected = correct_bgc_id_overview(mapping_file, mapping_path, genome_id)
        merged_dict.update(corrected)

    df = pd.DataFrame.from_dict(merged_dict).T
    df.index.name = "bgc_id"

    logging.debug(f"Writing file to: {table}")

    # Save dataframes to csv tables
    df_table = Path(table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(table)
    logging.info("Job done")
    return None


if __name__ == "__main__":
    gather_bgc_overview(sys.argv[1], sys.argv[2], sys.argv[3])
