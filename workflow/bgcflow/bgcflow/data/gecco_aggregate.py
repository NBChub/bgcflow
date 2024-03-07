import argparse
import logging
from pathlib import Path

import pandas as pd


def aggregate_gecco_cluster(file_list, outfile, cwd=".", sep="\t", suffix=".clusters"):
    """
    Aggregate GECCO data from multiple files into a single CSV file.

    This function reads a list of files from `file_list`, each file is expected to contain GECCO data.
    It then processes each file, adds a new column "genome_id" which is derived from the file name,
    and concatenates all data into a single pandas DataFrame. The resulting DataFrame is then saved
    to `outfile` in CSV format.

    Parameters:
    file_list (str): Path to a file containing a list of files to process. Each line in this file should be a path to a file.
    cwd (str): The current working directory. This is used to resolve relative paths in `file_list`.
    outfile (str): Path where the output CSV file will be saved.
    sep (str, optional): The column separator used in the input files. Defaults to "\t".
    suffix (str, optional): A suffix to strip from the file names when generating the "genome_id" column. Defaults to ".clusters".

    Returns:
    None
    """
    logging.info(f"Processing GECCO data from {file_list}...")

    file_list = Path(file_list)
    cwd = Path(cwd)
    gecco_summary_json = {
        (cwd / i).stem.strip(suffix): pd.read_csv((cwd / i), sep=sep)
        for i in pd.read_csv(file_list, header=None).loc[:, 0]
    }
    for k, v in gecco_summary_json.items():
        v["genome_id"] = k

    df = pd.concat(gecco_summary_json.values())

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(outfile, index=False)

    logging.info(f"GECCO data processed and saved to {outfile}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate GECCO data from multiple files into a single CSV file."
    )
    parser.add_argument(
        "file_list", help="Path to the file containing the list of files to process."
    )
    parser.add_argument("outfile", help="Path to the output file.")
    args = parser.parse_args()

    log_format = "%(levelname)-8s %(asctime)s   %(message)s"
    date_format = "%d/%m %H:%M:%S"
    logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

    aggregate_gecco_cluster(args.file_list, args.outfile)
