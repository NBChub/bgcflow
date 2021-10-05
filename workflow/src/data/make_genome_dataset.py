import os # Directory and file management
import pandas as pd # Dataframe
#import seaborn as sns # Visualization
#import matplotlib.pyplot as plt # Visualization 
#from matplotlib.ticker import StrMethodFormatter 
from Bio import SeqIO # Genbank processing 
from Bio.SeqUtils import GC

from pathlib import Path

def write_table(antismash_dir, genome_table, bgc_table):
    '''
    write df_genomes.csv and df_bgc_products.csv
    '''
    # Input 
    path = Path(antismash_dir)
    # Generate dataframe
    print('Generating dataframe for all genomes in the folder', path, '...')
    df_genomes, df_bgc_products = init_genome_dataframe(path)

    # Save dataframes to csv tables
    print('Saving dataframes to tables...')
    df_genomes.to_csv(genome_table)
    df_bgc_products.to_csv(bgc_table)
    return

def init_genome_dataframe(antismash_dir):
    '''
    Returns genome dataframe from input directory containing all antiSMASH results for the genomes
    '''
    path = Path(antismash_dir)
    genome_ids_list = [folder for folder in os.listdir(path) if os.path.isdir(path / folder)]
    column_names = ['genome_name', 'organism', 'closest_strain', 'ANI_closest_strain','genome_len', 'gc_content', 'accession_nbc', 
                    'topology', 'genus', 'species','records', 'plasmids', 'gene_count', 'bgcs_count',
                    'protoclusters_count', 'cand_clusters_count', 'date']
    
    df_genomes = pd.DataFrame(index=genome_ids_list, columns=column_names)
    df_bgc_products = pd.DataFrame(index=genome_ids_list)
    
    genomes_processed = 0
    for genome_id in genome_ids_list:
        gbk_file_path = os.path.join(antismash_dir, genome_id, genome_id + '.gbk')
        print('Genomes processed:', genomes_processed, 'Current genome:', genome_id)
        records_list = SeqIO.parse(gbk_file_path, 'genbank')
        
        record_cntr = 0
        gene_cntr = 0
        bgc_cntr = 0
        protoclusters_cntr = 0
        cand_clusters_cntr = 0
        bgc_type_dict = dict()
        A_count = 0
        C_count = 0
        G_count = 0
        T_count = 0
        length = 0
        
        for rec in records_list:
            record_cntr = record_cntr + 1
            # Information about organism 
            df_genomes.loc[genome_id, 'genome_name'] = rec.description
            df_genomes.loc[genome_id, 'organism'] = rec.annotations['organism']
            
            # Text mining to get closest organism and the ANI 
            # When source is of the format Streptomyces sp._root369_0.9634 # Need to fix this format
            source = rec.annotations['source']
            
            #if 'Unknown' in source:
            #    df_genomes.loc[genome_id, 'closest_strain'] = 'Unknown'
            #    df_genomes.loc[genome_id, 'ANI_closest_strain'] = 0
            #else:
            #    df_genomes.loc[genome_id, 'closest_strain'] = source.split('_0.')[0]
            #    df_genomes.loc[genome_id, 'ANI_closest_strain'] = float(source.split('_')[-1])
            
            df_genomes.loc[genome_id, 'topology'] = rec.annotations['topology']
            df_genomes.loc[genome_id, 'date'] = rec.annotations['date']
            df_genomes.loc[genome_id, 'accession_nbc'] = rec.id
            df_genomes.loc[genome_id, 'genus'] = rec.description.split(' ')[0]
            df_genomes.loc[genome_id, 'species'] = rec.description.split(' ')[1].split('_')[0]
            
            # Information on size and GC content
            A_count = A_count + rec.seq.count('A')
            C_count = C_count + rec.seq.count('C')
            G_count = G_count + rec.seq.count('G')
            T_count = T_count + rec.seq.count('T')
            length = length + len(rec.seq)
                        
            # Information on the number of BGCs, protoclusters and candidate clusters
            # Number of BGCs per type 
            # Number of genes
            for feat in rec.features:
                if feat.type == 'CDS':
                    gene_cntr = gene_cntr + 1
                if feat.type == 'region':
                    bgc_cntr = bgc_cntr + 1
                    bgc_type = '.'.join(sorted(feat.qualifiers['product']))
                    if bgc_type in bgc_type_dict.keys():
                        bgc_type_dict[bgc_type] = bgc_type_dict[bgc_type] + 1
                    else:
                        bgc_type_dict[bgc_type] = 1
                
                if feat.type == 'protocluster':
                    protoclusters_cntr = protoclusters_cntr + 1
                if feat.type == 'cand_cluster':
                    cand_clusters_cntr = cand_clusters_cntr + 1
                                    
        df_genomes.loc[genome_id, 'gene_count'] = gene_cntr
        df_genomes.loc[genome_id, 'bgcs_count'] = bgc_cntr
        df_genomes.loc[genome_id, 'protoclusters_count'] = protoclusters_cntr
        df_genomes.loc[genome_id, 'cand_clusters_count'] = cand_clusters_cntr
              
        for bgc_type in bgc_type_dict.keys():
            # Add column for BGC type if not present
            if bgc_type not in df_bgc_products.columns:
                df_bgc_products[bgc_type] = 0
            
            df_bgc_products.loc[genome_id, bgc_type] = bgc_type_dict[bgc_type]
        
        # Number of records in the genome
        #### Update plasmid counter: Currenly assumed records - 1
        df_genomes.loc[genome_id, 'records'] = record_cntr
        df_genomes.loc[genome_id, 'plasmids'] = record_cntr - 1 
        
        genomes_processed = genomes_processed + 1
                
        # Save genome len and GC content to dataframe
        df_genomes.loc[genome_id, 'gc_content'] = float(C_count + G_count) / length
        df_genomes.loc[genome_id, 'genome_len'] = length
    
    # Ordering df_bgc_products columns according to abundance 
    df_bgc_sum = df_bgc_products.sum(0)
    df_bgc_sum.sort_values(ascending=False, inplace=True)
    df_bgc_products = df_bgc_products.reindex(columns=df_bgc_sum.index)

    return df_genomes, df_bgc_products

#print(snakemake.input.antismash_dir, snakemake.output.df_genomes, snakemake.output.df_bgc_products)
write_table(snakemake.input.antismash_dir, snakemake.output.df_genomes, snakemake.output.df_bgc_products)