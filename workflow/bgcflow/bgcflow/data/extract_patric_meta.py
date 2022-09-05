import os, sys
import pandas as pd

def extract_ncbi_information(ncbi_meta_path, patric_genome_summary, patric_genome_metadata, patric_meta_path):
    """
    Extract PATRIC assembly json reports generated by get_assembly_information.py in a given directory into a csv file 

    Parameters
    ---------- 
    1. assembly_report_path : str / path
        Three different input types can be used: 
        - A directory containing assembly information json files generated by get_assembly_information.py
        - A single json file generated by get_assembly_information.py
        - A list of files as string separated by spaces (put the string inside '' in bash expression)
    2. outfile : str / path
        Location of the output csv file

    Returns
    -------
    1. outfile : .csv file
        A comma separated table summarizing the NCBI assembly reports metadata
    '''
    """

    df_ncbi = pd.read_csv(ncbi_meta_path, index_col="genome_id")
    df_patric_genome_summary = pd.read_csv(patric_genome_summary, sep="\t", index_col="genome_id")
    df_patric_genome_metadata = pd.read_csv(patric_genome_metadata, sep="\t", index_col="genome_id")

    df_patric_genome_metadata = df_patric_genome_metadata[df_patric_genome_metadata["genome_status"] != 'Plasmid']
    
    df_patric_meta = pd.DataFrame(index=df_ncbi.index)

    for genome_id in df_ncbi.index:
        refseq = df_ncbi.loc[genome_id, "refseq"]
        genbank = df_ncbi.loc[genome_id, "genbank"]
        if refseq in df_patric_genome_metadata["assembly_accession"].tolist():
            patric_genome_ids_list = df_patric_genome_metadata[df_patric_genome_metadata["assembly_accession"] == refseq].index.tolist()
            for patric_genome_id in patric_genome_ids_list:
                df_patric_meta.loc[genome_id, 'patric_genome_id'] = patric_genome_id
                for col in df_patric_genome_metadata.columns:
                    df_patric_meta.loc[genome_id, col] = df_patric_genome_metadata.loc[patric_genome_id, col]
                if patric_genome_id in df_patric_genome_summary.index:
                    for col in df_patric_genome_summary.columns:
                        df_patric_meta.loc[genome_id, col] = df_patric_genome_summary.loc[patric_genome_id, col]
        elif genbank in df_patric_genome_metadata["assembly_accession"].tolist():
            patric_genome_ids_list = df_patric_genome_metadata[df_patric_genome_metadata["assembly_accession"] == genbank].index.tolist()
            for patric_genome_id in patric_genome_ids_list:
                df_patric_meta.loc[genome_id, 'patric_genome_id'] = patric_genome_id
                df_patric_meta.loc[genome_id, :] = df_patric_genome_metadata.loc[patric_genome_id, :]
                if patric_genome_id in df_patric_genome_summary.index:
                    for col in df_patric_genome_summary.columns:
                        df_patric_meta.loc[genome_id, col] = df_patric_genome_summary.loc[patric_genome_id, col]
    df_patric_meta.index.name = 'genome_id'
    df_patric_meta.to_csv(patric_meta_path)

    return None

if __name__ == "__main__":
    extract_ncbi_information(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])