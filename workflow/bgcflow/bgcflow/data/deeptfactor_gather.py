import pandas as pd
from pathlib import Path
from alive_progress import alive_bar
import json
import sys
import logging

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def combine_deeptfactor_prediction(input_json, filter_str='_deeptfactor'):
    container = {}
    logging.info(f"Reading json files...")
    with alive_bar(len(input_json), title='Merging DeepTFactor prediction:') as bar:
        for item in input_json:
            item = Path(item)
            genome_id = item.stem
            if filter_str in genome_id:
                genome_id = genome_id.replace(filter_str, '')
            logging.debug(f"Processing {genome_id}")
            with open(item, "r") as f:
                data = json.load(f)
                container.update(data)
            bar()
    return container

def write_deeptf_table(input_json, deeptf_table):
    '''
    Write df_deeptfactor.csv table in processed data
    '''
    # Handle multiple json
    input_json = input_json.split()
    df = combine_deeptfactor_prediction(input_json)
    df = pd.DataFrame.from_dict(df).T
    df.index.name = "locus_tag"

    logging.debug(f"Writing file to: {deeptf_table}")

    # Save dataframes to csv tables
    df_table = Path(deeptf_table)
    df_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(deeptf_table, index=True)
    logging.info(f"Job done")
    return None

if __name__ == "__main__":
    write_deeptf_table(sys.argv[1], sys.argv[2])
