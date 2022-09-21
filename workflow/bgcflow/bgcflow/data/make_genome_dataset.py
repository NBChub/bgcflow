import os
import sys
import pandas as pd 
from Bio import SeqIO
from alive_progress import alive_bar
import logging
from pathlib import Path

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def write_genome_table(fna_dir, antismash_dir, samples_table, genome_table):
    '''
    Write df_genomes.csv table in processed data
    '''
    # Accomodate multiple inputs to generate dataframe
    shell_input = samples_table.split()
    logging.info(f"Reading samples table: {shell_input}")
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)

    # Generate dataframe
    df_genomes = update_bgc_info(antismash_dir, df_samples)

    # Save dataframes to csv tables
    genome_table = Path(genome_table)
    genome_table.parent.mkdir(parents=True, exist_ok=True)
    df_genomes.to_csv(genome_table)

    return None

def update_bgc_info(antismash_dir, df_genomes):
    '''
    Expands genome dataframe from input directory containing all antiSMASH results for the genomes
    '''
    # container for bgc products
    df_bgc_products = {}
    df_bgc_counter = {}

    # iterate over antismash genbanks
    logging.info(f"Parsing AntiSMASH results...")
    with alive_bar(df_genomes.shape[0], title='Updating BGC information:') as bar:
        for genome_id in df_genomes.index:
            logging.debug(f"Summarize BGCs info for {genome_id}")
            gbk_file_path = os.path.join(antismash_dir, genome_id, genome_id + '.gbk')
            records_list = SeqIO.parse(gbk_file_path, 'genbank')

            # statistics counter
            bgc_stats = {}
            bgc_cntr = 0
            protoclusters_cntr = 0
            cand_clusters_cntr = 0
            contig_edge_cntr = 0

            # product capture
            bgc_type_dict = {}

            for rec in records_list:
                # Information on the number of BGCs, protoclusters and candidate clusters
                for feat in rec.features:
                    if feat.type == 'region':
                        # get region counts
                        bgc_cntr = bgc_cntr + 1
                        if feat.qualifiers['contig_edge'][0] == 'True':
                            contig_edge_cntr = contig_edge_cntr + 1

                        # get product counts
                        bgc_type = '.'.join(sorted(feat.qualifiers['product']))
                        if bgc_type in bgc_type_dict.keys():
                            bgc_type_dict[bgc_type] = bgc_type_dict[bgc_type] + 1
                        else:
                            bgc_type_dict[bgc_type] = 1

                    if feat.type == 'protocluster':
                        protoclusters_cntr = protoclusters_cntr + 1
                    if feat.type == 'cand_cluster':
                        cand_clusters_cntr = cand_clusters_cntr + 1

            bgc_stats['bgcs_count'] = bgc_cntr
            bgc_stats['bgcs_on_contig_edge'] = contig_edge_cntr
            bgc_stats['protoclusters_count'] = protoclusters_cntr
            bgc_stats['cand_clusters_count'] = cand_clusters_cntr

            df_bgc_products[genome_id] = bgc_type_dict
            df_bgc_counter[genome_id] = bgc_stats
            bar()

        logging.debug(f"Merging information...")
        df_bgc_counter = pd.DataFrame.from_dict(df_bgc_counter).T
        df_bgc_products = pd.DataFrame.from_dict(df_bgc_products).T

        df_genomes = pd.concat([df_genomes, df_bgc_counter, df_bgc_products], axis=1)
        logging.debug(f"Done!")


    return df_genomes

if __name__ == "__main__":
    write_genome_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])