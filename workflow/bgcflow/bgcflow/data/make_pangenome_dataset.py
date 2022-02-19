import os
import sys
import pandas as pd
from shutil import copyfile

def get_roary_data(roary_interim_folder, roary_processed_folder):
    '''
    Copy important files from Roary interim to proecessed directory 

    Parameters
    ---------- 
    1. roary_folder : str / path
        Location of the output from Roary in the interim directory
    
    Returns
    -------
    1. roary_processed_folder : str / path
        Location of the processed output directory for important Roary results
    '''
    
    gene_presence_path = os.path.join(roary_interim_folder, 'gene_presence_absence.csv')
    df_gene_presence_summary = pd.read_csv(gene_presence_path, index_col='Gene')
    
    # Extract gene annotation columns to separate dataframe
    gene_summary_columns = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences',
        'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment',
        'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
        'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']

    gene_summary_out = os.path.join(roary_processed_folder, 'pangene_summary.csv')
    df_gene_summary = df_gene_presence_summary[gene_summary_columns].fillna('')
    df_gene_summary.to_csv(gene_summary_out)

    # Extract locus tags 
    df_gene_presence = df_gene_presence_summary.drop(columns=gene_summary_columns).fillna('')
    gene_presence_out = os.path.join(roary_processed_folder, 'gene_presence_absence.csv')
    df_gene_presence.to_csv(gene_presence_out)

    # Save the gene presence absence binary matrix
    gene_presence_binary_path = os.path.join(roary_interim_folder, 'gene_presence_absence.Rtab')
    gene_presence_binary_out = os.path.join(roary_processed_folder, 'gene_presence_absence_binary.csv')
    df_gene_presence_binary = pd.read_csv(gene_presence_binary_path, index_col='Gene', sep='\t')
    df_gene_presence_binary.to_csv(gene_presence_binary_out)

    # Copy other output files to processed directory
    copyfile(os.path.join(roary_interim_folder, 'conserved_vs_total_genes.png'),
        os.path.join(roary_processed_folder, 'conserved_vs_total_genes.png'))

    copyfile(os.path.join(roary_interim_folder, 'Rplots.pdf'),
        os.path.join(roary_processed_folder, 'Rplots.pdf'))

    copyfile(os.path.join(roary_interim_folder, 'pan_genome_reference.fa'),
        os.path.join(roary_processed_folder, 'pan_genome_reference.fa'))
    
    copyfile(os.path.join(roary_interim_folder, 'summary_statistics.txt'),
        os.path.join(roary_processed_folder, 'summary_statistics.txt'))

    copyfile(os.path.join(roary_interim_folder, 'number_of_conserved_genes.Rtab'),
        os.path.join(roary_processed_folder, 'number_of_conserved_genes.Rtab'))
    
    copyfile(os.path.join(roary_interim_folder, 'number_of_genes_in_pan_genome.Rtab'),
        os.path.join(roary_processed_folder, 'number_of_genes_in_pan_genome.Rtab'))
    
    return df_gene_presence

if __name__ == "__main__":
    get_roary_data(sys.argv[1], sys.argv[2])