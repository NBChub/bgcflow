import os
import pandas as pd

from pathlib import Path

def write_ncbi_meta(assembly_report_path, meta_out_path):
    '''
    Write df_ncbi_meta.csv with metadata from NCBI assembly reports
    '''

    # Input 
    path = Path(assembly_report_path)
    # Generate dataframe
    print('Generating dataframe for NCBI assembly metadata from folder', path, '...')
    df_ncbi_meta = get_assembly_meta(path)

    # Save dataframes to csv tables
    print('Saving metadata dataframes to tables...')
    df_ncbi_meta.to_csv(meta_out_path)

    return None


def get_assembly_meta(assembly_report_path):
    '''
    Returns metadata dataframe from input directory containing all assembly reports downloaded using 
    ncbi-genome-download https://github.com/kblin/ncbi-genome-download
    Use assembly_report as --formats while running ncbi-genome-download 
    '''
    
    path = Path(assembly_report_path)
    
    genome_list = [file[:-4] for file in os.listdir(path) if '.txt' in file]
    
    # List of columns in df_ncbi_meta
    ncbi_meta_columns = ['assembly', 'organism', 'genus', 'species', 'strain', 'tax_id', 'refseq_category', 'refseq', 'genbank', 'refseq_genbank_identity', 'biosample', 'submitter', 'date']
    
    df_ncbi_meta = pd.DataFrame(index = genome_list, columns = ncbi_meta_columns)
    df_ncbi_meta.index.name = 'genome_id'
    
    for genome_id in genome_list:
        report_path = os.path.join(path, genome_id + '.txt')
        df_meta_genome = pd.read_csv(report_path, sep=':', header=None, index_col=0)
        assembly_accn = df_meta_genome.loc['# Assembly name', 1].strip()
            
        if assembly_accn.startswith('ASM'):
            df_ncbi_meta.loc[genome_id, 'assembly'] = assembly_accn
        else:
            print('Assembly accession not in standard format', assembly_accn, genome_id)
        df_ncbi_meta.loc[genome_id, 'assembly'] = assembly_accn
        df_ncbi_meta.loc[genome_id, 'organism'] = df_meta_genome.loc['# Organism name', 1].strip()
        df_ncbi_meta.loc[genome_id, 'genus'] = df_meta_genome.loc['# Organism name', 1].strip().split(' ')[0]
        df_ncbi_meta.loc[genome_id, 'species'] = df_meta_genome.loc['# Organism name', 1].strip().split(' ')[1]
        if '# Infraspecific name' in df_meta_genome.index:
            strain = df_meta_genome.loc['# Infraspecific name', 1].strip()
            if 'strain=' in strain:
                df_ncbi_meta.loc[genome_id, 'strain'] = strain.split('=')[1]
            elif '# Isolate' in df_meta_genome.index:
                strain = df_meta_genome.loc['# Isolate', 1].strip()
                df_ncbi_meta.loc[genome_id, 'strain'] = strain
                
        df_ncbi_meta.loc[genome_id, 'tax_id'] = df_meta_genome.loc['# Taxid', 1].strip()
        df_ncbi_meta.loc[genome_id, 'biosample'] = df_meta_genome.loc['# BioSample', 1].strip()
        df_ncbi_meta.loc[genome_id, 'submitter'] = df_meta_genome.loc['# Submitter', 1].strip()
        df_ncbi_meta.loc[genome_id, 'date'] = df_meta_genome.loc['# Date', 1].strip()
            
        if '# RefSeq category' in df_meta_genome.index:
            df_ncbi_meta.loc[genome_id, 'refseq_category'] = df_meta_genome.loc['# RefSeq category', 1].strip()
        if '# RefSeq assembly accession' in df_meta_genome.index:
            df_ncbi_meta.loc[genome_id, 'refseq'] = df_meta_genome.loc['# RefSeq assembly accession', 1].strip()
        if '# GenBank assembly accession' in df_meta_genome.index:
            df_ncbi_meta.loc[genome_id, 'genbank'] = df_meta_genome.loc['# GenBank assembly accession', 1].strip()
        if '# RefSeq assembly and GenBank assemblies identical' in df_meta_genome.index:
            df_ncbi_meta.loc[genome_id, 'refseq_genbank_identity'] = df_meta_genome.loc['# RefSeq assembly and GenBank assemblies identical', 1].strip()
            
    return df_ncbi_meta

write_ncbi_meta(snakemake.input.assembly_report_path, snakemake.output.meta_out_path)