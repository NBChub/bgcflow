from pathlib import Path
import sys
import pandas as pd
from alive_progress import alive_bar
import json
import logging

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def gather_seqfu_json(seqfu_input, output_file):
    """
    Merge seqfu json input into one csv table
    """
    # Accomodate multiple inputs to generate dataframe
    shell_input = seqfu_input.split()
    logging.info(f"Merging seqfu results...")
    container = {}
    with alive_bar(len(shell_input), title="Updating BGC information:") as bar:
        for d in shell_input:
            logging.debug(f"Reading input: {Path(d).name}")
            with open(d, "r") as f:
                read = json.load(f)
                genome_id = read[0]["Filename"]
                read[0].pop("Filename")
                container[genome_id] = read[0]
            bar()
    logging.info(f"Converting dictionary to table...")
    df = pd.DataFrame.from_dict(container).T
    df.index.name = "genome_id"
    df.to_csv(output_file)
    logging.info(f"Job done!")
    return


if __name__ == "__main__":
    gather_seqfu_json(sys.argv[1], sys.argv[2])
