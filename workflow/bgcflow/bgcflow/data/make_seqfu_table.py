import json
import logging
import sys
from pathlib import Path

import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def gather_seqfu_json(seqfu_input, output_file):
    """
    Merge seqfu json input into one csv table
    """
    # Accomodate multiple inputs to generate dataframe
    shell_input = Path(seqfu_input)
    logging.info(shell_input)
    if shell_input.is_file() and shell_input.suffix == ".json":
        logging.info(f"Getting SeqFu overview from a single file: {shell_input}")
        shell_input_files = shell_input

    elif shell_input.is_file() and shell_input.suffix == ".txt":
        logging.info(f"Getting SeqFu overview  from a text file: {shell_input}")
        with open(shell_input, "r") as file:
            file_content = [i.strip("\n") for i in file.readlines()]
            if len(file_content) == 1:
                # Paths space-separated on a single line
                paths = file_content[0].split()
            else:
                # Paths written on separate lines
                paths = file_content
            shell_input_files = [
                Path(path) for path in paths if Path(path).suffix == ".json"
            ]
    else:
        shell_input_files = [
            Path(file)
            for file in str(shell_input).split()
            if Path(file).suffix == ".json"
        ]
        logging.info(
            f"Getting SeqFu overview from the given list of {len(shell_input_files)} files..."
        )

    logging.info("Merging SeqFu results...")
    container = {}
    with alive_bar(len(shell_input_files), title="Updating SeqFu information:") as bar:
        for d in shell_input_files:
            logging.debug(f"Reading input: {Path(d).name}")
            with open(d, "r") as f:
                read = json.load(f)
                genome_id = read[0]["Filename"]
                read[0].pop("Filename")
                container[genome_id] = read[0]
            bar()
    logging.info("Converting dictionary to table...")
    df = pd.DataFrame.from_dict(container).T
    df.index.name = "genome_id"
    df.to_csv(output_file)
    logging.info("Job done!")
    return


if __name__ == "__main__":
    gather_seqfu_json(sys.argv[1], sys.argv[2])
