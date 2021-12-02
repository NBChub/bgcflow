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

        with open(report_path, 'r') as report_file:
            lines = report_file.readlines()
            strain_found = False

            for line in lines:
                if line.startswith('# Assembly name:'):
                    assembly_accn = line.split('Assembly name:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'assembly'] = assembly_accn
                elif line.startswith('# Organism name'):
                    organism = line.split('Organism name:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'organism'] = organism
                    df_ncbi_meta.loc[genome_id, 'genus'] = organism.strip().split(' ')[0]
                    df_ncbi_meta.loc[genome_id, 'species'] = organism.strip().split(' ')[1]
                elif line.startswith('# Infraspecific name:'):
                    strain = line.split('Infraspecific name:')[1].strip()
                    if 'strain=' in strain:
                        df_ncbi_meta.loc[genome_id, 'strain'] = strain.split('strain=')[1]
                        strain_found = True
                elif line.startswith('# Isolate:'):
                    if strain_found == False:
                        strain = line.split('Isolate:')[1].strip()
                        df_ncbi_meta.loc[genome_id, 'strain'] = strain
                        strain_found = True
                elif line.startswith('# Taxid'):
                    tax_id = line.split('Taxid:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'tax_id'] = tax_id
                elif line.startswith('# BioSample'):
                    biosample = line.split('BioSample:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'biosample'] = biosample
                elif line.startswith('# Submitter'):
                    submitter = line.split('Submitter:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'submitter'] = submitter
                elif line.startswith('# Date'):
                    date = line.split('Date:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'date'] = date
                elif line.startswith('# RefSeq category'):
                    refseq_category = line.split('RefSeq category:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'refseq_category'] = refseq_category
                elif line.startswith('# RefSeq assembly accession'):
                    refseq = line.split('RefSeq assembly accession:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'refseq'] = refseq
                elif line.startswith('# GenBank assembly accession'):
                    genbank = line.split('GenBank assembly accession:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'genbank'] = genbank
                elif line.startswith('# RefSeq assembly and GenBank assemblies identical'):
                    refseq_genbank_identity = line.split('RefSeq assembly and GenBank assemblies identical:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'refseq_genbank_identity'] = refseq_genbank_identity
                
            if strain_found == False:
                df_ncbi_meta.loc[genome_id, 'strain'] = genome_id

    return df_ncbi_meta

write_ncbi_meta(snakemake.input.assembly_report_path, snakemake.output.meta_out_path)