import argparse
import pandas as pd
import os

def remove_special_char(s):
    if '/' in str(s):
        return(s.replace('/','_'))
    elif "'" in str(s):
        return(s.replace("'",'_variant'))
    elif "(" in str(s):
        s = s.replace("(",'_')
        return(s.replace(")",''))
    else:
        return(s)

def process_genus_species(genus, species):
    """
    Process the given genus and species, read the corresponding CSV files, apply the remove_special_char function
    to the 'Gene' column of each DataFrame, and then save the modified DataFrames back to their original CSV files.

    Parameters:
    genus (str): The name of the genus.
    species (str): The name of the species.

    Returns:
    None
    """
    print(f'Start getting summary table:  {genus}, {species}')

    apm_binary = pd.read_csv('../data/genus/' + genus +'/'+ species + '/roary/df_gene_presence_binary.csv', low_memory = False)
    summary = pd.read_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_summary.csv',  low_memory = False)
    summary_2 = pd.read_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_summary_v2.csv',  low_memory = False)
    annotation = pd.read_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_eggnog_summary.csv',  low_memory = False)
    locustag = pd.read_csv('../data/genus/' + genus +'/'+ species + '/roary/df_gene_presence_locustag.csv',  low_memory = False)

    apm_binary.Gene = apm_binary.Gene.apply(remove_special_char)
    summary.Gene = summary.Gene.apply(remove_special_char)
    summary_2.Gene = summary_2.Gene.apply(remove_special_char)
    annotation.Gene = annotation.Gene.apply(remove_special_char)
    locustag.Gene = locustag.Gene.apply(remove_special_char)

    apm_binary.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_gene_presence_binary.csv', index = False)
    summary.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_summary.csv',  index = False)
    summary_2.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_summary_v2.csv',  index = False)
    annotation.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_eggnog_summary.csv',  index = False)
    locustag.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_gene_presence_locustag.csv',  index = False)

    print('Finished')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process genus name.')
    parser.add_argument('genus', type=str, help='Name of the genus')
    parser.add_argument('species', type=str, help='Name of the species')
    args = parser.parse_args()

    process_genus_species(args.genus, args.species)
