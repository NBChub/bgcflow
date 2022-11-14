import json
import pandas as pd
from pathlib import Path
import numpy as np
import sys

import logging
log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.INFO)

# regions
def region_table_builder(f, accession):
    """
    Given a feature of a record, return value to build database table
    """
    # grab region values
    region_number = f['qualifiers']['region_number'][0]
    region_id = f"{accession}.region{str(region_number).zfill(3)}"
    location = f['location'].strip("[").strip("]")
    start_pos, end_pos = location.split(":")
    contig_edge = f['qualifiers']['contig_edge'][0]
    # fill values
    value = {region_id : {'accession' : accession,
                          'region_number' : region_number,
                          'location' : location,
                          'start_pos' : start_pos,
                          'end_pos' : end_pos,
                          'contig_edge' : contig_edge,
                          'product' : f['qualifiers']['product'],
                          'rules' : f['qualifiers']['rules']
                         }
            }
    return value

def cdss_table_builder(f, cds_id):
    #location, strand = f['location'].strip("[").strip(")").split("](")
    try:
        gene_function = f['qualifiers']['gene_functions']
    except KeyError:
        logging.debug(f"{cds_id} does not have gene_function. Available values:\n{f['qualifiers'].keys()}")
        gene_function = None
        
    try:
        name = f['qualifiers']['gene'][0]
    except KeyError:
        logging.debug(f"{cds_id} does not have gene name. Available values:\n{f['qualifiers'].keys()}")
        name = None
    
    try:
        gene_kind = f['qualifiers']['gene_kind'][0]
    except KeyError:
        logging.debug(f"{cds_id} does not have gene kind. Available values:\n{f['qualifiers'].keys()}")
        gene_kind = None
        
    try:
        EC_number = f['qualifiers']['EC_number'][0]
    except KeyError:
        logging.debug(f"{cds_id} does not have EC_number. Available values:\n{f['qualifiers'].keys()}")
        EC_number = None
    
    value = {cds_id : {'gene_function' : gene_function,
                       'locus_tag' : f['qualifiers']['locus_tag'][0],
                       'name' : name,
                       'product' : f['qualifiers']['product'][0],
                       'translation' : f['qualifiers']['translation'][0],
                       'location' : f['location'],
                       'gene_kind' : gene_kind,
                       'codon_start' : f['qualifiers']['codon_start'][0],
                       'EC_number' : EC_number
                      }
            }
    return value

def region_finder(cdss_id, location, regions_container):
    location, strand = location.strip("[").strip(")").split("](")
    q_start, q_stop = [int(i) for i in location.split(":")]
    #print(cdss_id, q_start, q_stop, strand)
    query = range(q_start, q_stop)
    
    hits = []
    
    for region_id in regions_container.keys():
        start_pos = int(regions_container[region_id]['start_pos'])
        end_pos = int(regions_container[region_id]['end_pos'])
        target = range(start_pos, end_pos)
        try:
            value = range(max(query[0], target[0]), min(query[-1], target[-1])+1)
            if len(value) > 1:
                logging.debug(f"{cdss_id} overlaps with {region_id} at {value}")
                hits.append(region_id)
            #else:
            #    hits.append(np.nan)
        except IndexError:
            pass
    return hits

def antismash_json_exporter(json_file, output_dir):

    file = Path(json_file)

    table_dna_sequences = {}
    table_regions = {}
    table_cdss = {}

    with open(str(file), "r") as f:
        gbk = json.load(f)

    genome_id = Path(gbk['input_file']).stem
    logging.info(f"Extracting information from {gbk['input_file']}")
    for record in gbk['records']:
        sequence_id = record['id']
        logging.info(f"Getting dna_sequences information for {sequence_id}")

        # dna_sequences
        record_container = {}

        ## print(record.keys())
        ## print(2, record['seq']) # dna
        record_container['seq'] = record['seq']['data']

        ## print(3, record['description']) # definition
        record_container['description'] = record['description']

        ## print(5, {'molecule_type' : record['annotations']['molecule_type']}) # contig type?
        record_container['molecule_type'] = record['annotations']['molecule_type']

        ## print(5, {'topology' : record['annotations']['topology']}) # chromosome type?
        record_container['topology'] = record['annotations']['topology']

        ## print(5, {'accessions' : record['annotations']['accessions']}) # acc
        assert len(record['annotations']['accessions']) == 1
        record_container['accessions'] = record['annotations']['accessions'][0]

        record_container['genome_id'] = genome_id

        ## print(1, record['id']) # sequence_id
        table_dna_sequences[sequence_id ] = record_container

        # ---------------------------------------------------
        # GRAB REGIONS AND CDSS TABLE
        logging.info(f"Extracting regions and cds information from {sequence_id}")
        cds_ctr = 1
        region_ctr = 1
        accession = record_container['accessions']
        for f in record['features']:
                # Fill region table
                if f['type'] == 'region':
                    logging.info(f"Processing region: {region_ctr}")
                    table_regions.update(region_table_builder(f, accession))
                    region_ctr = region_ctr + 1

                # fill cdss table
                if f['type'] == 'CDS':
                    cds_id = f"{accession}-cds_{cds_ctr}"
                    cdss_data = cdss_table_builder(f, cds_id)
                    region_hits = region_finder(cds_id, cdss_data[cds_id]["location"], table_regions)
                    if len(region_hits) > 0:
                        cdss_data[cds_id]['region_id'] = region_hits[0]
                    else:
                        cdss_data[cds_id]['region_id'] = None
                    table_cdss.update(cdss_data)
                    cds_ctr = cds_ctr + 1
        # ---------------------------------------------------

        #print(7, record['areas']) # antismash regions
        #print(8, record['modules'].keys()) # dna
    
    # write json files
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    target_jsons = {'dna_sequences' : table_dna_sequences, 
                    'regions' : table_regions, 
                    'cdss' : table_cdss}
    for k in target_jsons.keys():
        with open(output_dir / f"{genome_id}_{k}.json", "w") as output_file:
            json.dump(target_jsons[k], output_file, indent=2)
    return

if __name__ == "__main__":
    antismash_json_exporter(sys.argv[1], sys.argv[2])