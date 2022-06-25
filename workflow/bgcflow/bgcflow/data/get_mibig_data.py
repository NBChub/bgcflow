import os
from pathlib import Path
import sys
import pandas as pd
import logging
import json

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def extract_mibig_info(mibig_json_path, mibig_bgc_table):
    '''
    Returns mibig BGC information critical for the known bgc table from JSON files

    Parameters
    ----------
    1. mibig_json_path : str / path 
        Location of the resources folder with all MIBIG JSON files extracted

    Returns
    -------
    1. mibig_bgc_table: str / path
        Dataframe with MIBIG BGCs information
    2. mibig_compound_table: str / path
        Dataframe with compounds information found at MIBIG
    '''
    
    if not os.path.isdir(mibig_json_path):
        raise FileNotFoundError(f'No such file or directory: {mibig_json_path}')

    df_mibig_bgcs = pd.DataFrame(columns=['biosyn_class', 'compounds'])
    df_mibig_compounds = pd.DataFrame(columns=[])
    
    mibig_list = [mibig_file.split('.')[0] for mibig_file in os.listdir(mibig_json_path) if '.json' in mibig_file]
    for mibig_id in mibig_list:
        mibig_id_file = os.path.join(mibig_json_path, mibig_id + '.json')
        with open(mibig_id_file, 'r') as json_obj:
            mibig_data = json.load(json_obj)
            df_mibig_bgcs.loc[mibig_id, 'biosyn_class'] = ';'.join(mibig_data.get('cluster').get('biosyn_class')) 
            
            compounds_list = mibig_data.get('cluster').get('compounds')
            df_mibig_bgcs.loc[mibig_id, 'compounds'] = ';'.join([compound.get('compound') for compound in compounds_list]) 

            chem_acts_list = []
            for compound in mibig_data.get('cluster').get('compounds'):
                if 'chem_acts' in compound.keys():
                    chem_acts = compound.get('chem_acts') 
                    for activity in chem_acts:
                        if activity not in chem_acts_list:
                            chem_acts_list.append(activity)
            df_mibig_bgcs.loc[mibig_id, 'chem_acts'] = ';'.join(chem_acts_list)
            
            if 'accession' in mibig_data.get('cluster').get('loci').keys():
                df_mibig_bgcs.loc[mibig_id, 'accession'] = mibig_data.get('cluster').get('loci').get('accession')
            if 'completeness' in mibig_data.get('cluster').get('loci').keys():
                df_mibig_bgcs.loc[mibig_id, 'completeness'] = mibig_data.get('cluster').get('loci').get('completeness')
            if 'evidence' in mibig_data.get('cluster').get('loci').keys():
                df_mibig_bgcs.loc[mibig_id, 'evidence'] = ';'.join(mibig_data.get('cluster').get('loci').get('evidence'))
            if 'organism_name' in mibig_data.get('cluster').keys():
                df_mibig_bgcs.loc[mibig_id, 'organism_name'] = mibig_data.get('cluster').get('organism_name')
            if 'ncbi_tax_id' in mibig_data.get('cluster').keys():
                df_mibig_bgcs.loc[mibig_id, 'ncbi_tax_id'] = mibig_data.get('cluster').get('ncbi_tax_id')
            if 'publications' in mibig_data.get('cluster').keys():
                df_mibig_bgcs.loc[mibig_id, 'publications'] = ';'.join(mibig_data.get('cluster').get('publications'))
            
    df_mibig_bgcs.sort_index(inplace=True)
    df_mibig_bgcs.index.name = 'mibig_id'

    df_mibig_bgcs.to_csv(mibig_bgc_table)

    return df_mibig_bgcs

if __name__ == "__main__":
    extract_mibig_info(sys.argv[1], sys.argv[2])