import os
import sys
import pandas as pd

"""def merge_multiple_inputs(samples_path):
    if type(assembly_report_path) == list:
        dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in assembly_report_path]
        df = pd.concat(dfList, axis=0)
    return"""

def get_ncbi_meta(assembly_report_path, df_samples):
    '''
    Returns metadata dataframe from input directory containing all assembly reports downloaded using 
    ncbi-genome-download https://github.com/kblin/ncbi-genome-download
    Use assembly_report as --formats while running ncbi-genome-download 
    '''

    # Extract NCBI genomes from the df_samples
    genome_list = df_samples[df_samples.source.eq("ncbi")].genome_id.to_list()
    
    # List of columns in df_ncbi_meta
    ncbi_meta_columns = ['assembly', 'organism', 'genus', 'species', 'strain', 'tax_id', 
                            'refseq_category', 'refseq', 'genbank', 'assembly_type', 'release_type',
                            'assembly_level', 'genome_representation', 'refseq_genbank_identity', 
                            'biosample', 'submitter', 'date']
    
    df_ncbi_meta = pd.DataFrame(index = genome_list, columns = ncbi_meta_columns)
    df_ncbi_meta.index.name = 'genome_id'
    
    for genome_id in genome_list:
        report_path = os.path.join(assembly_report_path, genome_id + '.txt')

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
                elif line.startswith('# Assembly type:'):
                    assembly_type = line.split('Assembly type:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'assembly_type'] = assembly_type
                elif line.startswith('# Release type:'):
                    release_type = line.split('Release type:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'release_type'] = release_type
                elif line.startswith('# Assembly level:'):
                    assembly_level = line.split('Assembly level:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'assembly_level'] = assembly_level
                elif line.startswith('# Genome representation:'):
                    genome_representation = line.split('Genome representation:')[1].strip()
                    df_ncbi_meta.loc[genome_id, 'genome_representation'] = genome_representation
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
    

def extract_ncbi_org_info(prokka_dir, df_ncbi_meta):
    """
    Returns organism_info.txt with genus, species, strian info for prokka run inputs
    This function returns values from NCBI assembly reports
    """

    for idx in df_ncbi_meta.index:
        GENUS = df_ncbi_meta.loc[idx, 'genus']
        SPECIES = df_ncbi_meta.loc[idx, 'species']
        STRAIN_ID = df_ncbi_meta.loc[idx, 'strain']

        if not os.path.isdir(os.path.join(prokka_dir, idx)):
            os.mkdir(os.path.join(prokka_dir, idx))
        org_info_path = os.path.join(prokka_dir, idx, 'organism_info.txt')
        with open(org_info_path, 'w') as file_obj:
            file_obj.write(','.join([GENUS,  SPECIES, STRAIN_ID]))
    
    return None


def extract_samples_org_info(prokka_dir, df_samples):
    """
    Returns organism_info.txt with genus, species, strian info for prokka run inputs
    This function returns values from provided sample.csv excluding NCBI
    """

    genome_list = df_samples[~df_samples.source.eq("ncbi")].genome_id.to_list()

    for idx in genome_list:
        GENUS = df_samples.loc[idx, 'genus']
        SPECIES = df_samples.loc[idx, 'species']
        STRAIN_ID = df_samples.loc[idx, 'strain']

        if not os.path.isdir(os.path.join(prokka_dir, idx)):
            os.mkdir(os.path.join(prokka_dir, idx))
        org_info_path = os.path.join(prokka_dir, idx, 'organism_info.txt')
        with open(org_info_path, 'w') as file_obj:
            file_obj.write(','.join([GENUS,  SPECIES, STRAIN_ID]))
    
    return None


def extract_org_info(samples_path, assembly_report_path, prokka_dir, ncbi_meta_path):
    '''
    Write df_ncbi_meta.csv with metadata from NCBI assembly reports
    '''
    # wrap single or multiple inputs & generate dataframe
    shell_input = samples_path.split()
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)
    
    df_ncbi_meta = get_ncbi_meta(assembly_report_path, df_samples)

    # Save dataframes to csv tables
    df_ncbi_meta.to_csv(ncbi_meta_path)

    # Extract organism infor for prokka input
    extract_samples_org_info(prokka_dir, df_samples)
    extract_ncbi_org_info(prokka_dir, df_ncbi_meta)

    return None


if __name__ == "__main__":
    extract_org_info(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])