from pathlib import Path
from Bio import SeqIO
import sys
import pandas as pd
from alive_progress import alive_bar
import json
import logging

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def combine_bgc_counts(input_json, filter_str='_bgc_counts'):
    container = {}
    logging.info(f"Reading json files...")
    with alive_bar(len(input_json), title='Merging BGC counts:') as bar:
        for item in input_json:
            item = Path(item)
            genome_id = item.stem
            if filter_str in genome_id:
                genome_id = genome_id.replace(filter_str, '')
            logging.debug(f"Processing {genome_id}")
            with open(item, "r") as f:
                data = json.load(f)
                products = data[genome_id].pop("products")
                data[genome_id] = data[genome_id] | products
                container.update(data)
            bar()
    return container

def write_genome_table(input_json, samples_table, genome_table):
    '''
    Write df_genomes.csv table in processed data
    '''
    # Accomodate multiple inputs to generate dataframe
    shell_input = samples_table.split()
    logging.info(f"Reading samples table: {shell_input}")
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)

    # Handle multiple json
    input_json = input_json.split()
    bgc_counts = combine_bgc_counts(input_json)
    bgc_counts = pd.DataFrame.from_dict(bgc_counts).T

    logging.debug(f"Writing file to: {genome_table}")

    # Generate dataframe
    df_genomes = pd.concat([df_samples, bgc_counts], axis=1)

    # Save dataframes to csv tables
    genome_table = Path(genome_table)
    genome_table.parent.mkdir(parents=True, exist_ok=True)
    df_genomes.to_csv(genome_table, index=False)
    logging.info(f"Job done")
    return None

if __name__ == "__main__":
    write_genome_table(sys.argv[1], sys.argv[2], sys.argv[3])