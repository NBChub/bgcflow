import os
import sys
import pandas as pd 
from Bio import SeqIO
from alive_progress import alive_bar

from pathlib import Path

def write_genome_table(fna_dir, antismash_dir, samples_table, genome_table):
    '''
    Write df_genomes.csv table in processed data
    '''
    # Accomodate multiple inputs to generate dataframe
    shell_input = samples_table.split()
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)
    
    # Generate dataframe
    df_genomes = init_genome_dataframe(fna_dir, df_samples)
    df_genomes = update_bgc_info(antismash_dir, df_genomes)
    # Save dataframes to csv tables
    df_genomes.to_csv(genome_table)

    return None

def init_genome_dataframe(fna_dir, df_samples):
    '''
    Returns genome dataframe with GC content and number of contigs per genome
    '''

    df_genomes = df_samples.copy()
    df_genomes[['contigs', 'gc_content', 'genome_len']] = 0

    with alive_bar(df_genomes.shape[0]) as bar:
        for genome_id in df_genomes.index:
            fna_file_path = os.path.join(fna_dir, str(genome_id) + '.fna')
            records_list = SeqIO.parse(fna_file_path, 'fasta')
            
            record_cntr = 0
            A_count = 0
            C_count = 0
            G_count = 0
            T_count = 0
            length = 0
            
            for rec in records_list:
                record_cntr = record_cntr + 1
                # Information on size and GC content
                A_count = A_count + rec.seq.count('A')
                C_count = C_count + rec.seq.count('C')
                G_count = G_count + rec.seq.count('G')
                T_count = T_count + rec.seq.count('T')

                if 'A' not in rec.seq:
                    A_count = A_count + rec.seq.count('a')
                    C_count = C_count + rec.seq.count('c')
                    G_count = G_count + rec.seq.count('g')
                    T_count = T_count + rec.seq.count('t')
                length = length + len(rec.seq)            
                
            # Number of records in the genome
            df_genomes.loc[genome_id, 'contigs'] = int(record_cntr)
                            
            # Save genome len and GC content to dataframe
            df_genomes.loc[genome_id, 'gc_content'] = float(C_count + G_count) / length
            df_genomes.loc[genome_id, 'genome_len'] = int(length)

            bar()
            
        df_genomes = df_genomes.astype({'genome_len': 'int', 'contigs': 'int'})
    
    return df_genomes

def update_bgc_info(antismash_dir, df_genomes):
    '''
    Expands genome dataframe from input directory containing all antiSMASH results for the genomes
    '''

    for genome_id in df_genomes.index:
        gbk_file_path = os.path.join(antismash_dir, genome_id, genome_id + '.gbk')
        records_list = SeqIO.parse(gbk_file_path, 'genbank')

        bgc_cntr = 0
        protoclusters_cntr = 0
        cand_clusters_cntr = 0
        contig_edge_cntr = 0
        
        for rec in records_list:
            # Information on the number of BGCs, protoclusters and candidate clusters
            for feat in rec.features:
                if feat.type == 'region':
                    bgc_cntr = bgc_cntr + 1
                    if feat.qualifiers['contig_edge'] == 'True':
                        contig_edge_cntr = contig_edge_cntr + 1
                if feat.type == 'protocluster':
                    protoclusters_cntr = protoclusters_cntr + 1
                if feat.type == 'cand_cluster':
                    cand_clusters_cntr = cand_clusters_cntr + 1
                                    
        df_genomes.loc[genome_id, 'bgcs_count'] = bgc_cntr
        df_genomes.loc[genome_id, 'bgcs_on_contig_edge'] = contig_edge_cntr
        df_genomes.loc[genome_id, 'protoclusters_count'] = protoclusters_cntr
        df_genomes.loc[genome_id, 'cand_clusters_count'] = cand_clusters_cntr
              
    return df_genomes

if __name__ == "__main__":
    write_genome_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])