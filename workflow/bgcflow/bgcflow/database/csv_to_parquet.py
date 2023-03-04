import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.INFO)


def csv_to_parquet(project_folder, output_folder="."):
    """
    Given a list of csv files, convert into parquet files
    """
    project_folder = Path(project_folder)
    if output_folder == ".":
        output_folder = project_folder / "data_warehouse"
    else:
        pass

    logging.info(
        f"Grabbing all csv from folder {project_folder} and saving parquets in {output_folder}"
    )
    for i in project_folder.rglob("*.csv"):
        if "docs" in str(i):
            pass
        elif "dbt" in str(i):
            pass
        elif "ipynb_checkpoints" in str(i):
            pass
        else:
            category = str(i).split("/")[3]
            output_subfolder = output_folder / category
            df = pd.read_csv(i)
            output_parquet = output_subfolder / f"{i.stem}.parquet"
            logging.info(f"Converting {i} to {output_parquet}")
            output_subfolder.mkdir(parents=True, exist_ok=True)
            df.to_parquet(output_parquet)

    return


if __name__ == "__main__":
    csv_to_parquet(sys.argv[1])
