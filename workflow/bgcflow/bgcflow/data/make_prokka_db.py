import ncbi_genome_download as ngd
import pandas as pd
import os, sys, glob, gzip, io

def build_prokka_refgbff(name, prokka_db_table, outfile):
    """
    Given a project name and csv table path, download all NCBI accession id
    and returns a concatenated gbk file for prokka --proteins params.
    """
    download_dir = f"resources/prokka_db/{name}"
    #outfile = f"resources/prokka_db/reference_{name}.gbff"
    
    # Generate download directory
    os.makedirs(download_dir, exist_ok=True)
    
    # Get lists of accession ids and its sources
    df = pd.read_csv(prokka_db_table)
    ngd_input = {"refseq" : [],
                 "genbank" : []
                }
    for acc in df.Accession:
        if acc.startswith("GCF"):
            ngd_input['refseq'].append(acc)
        elif acc.startswith("GCA"):
            ngd_input['genbank'].append(acc)
        else:
            raise
    
    # Download gbff with NCBI genome download     
    for s in ngd_input.keys():
        if ngd_input[s]:
            acc_list = ",".join(ngd_input[s])
            ngd.download(section=s,
                         file_formats='genbank', 
                         assembly_accessions=acc_list,  
                         output=download_dir, 
                         groups="bacteria")
    
    # Concatenate gbff
    reference_files = glob.glob(f"{download_dir}/*/*/*/*.gbff.gz")

    with open(outfile, 'w') as outfile:
        for names in reference_files:
            with io.TextIOWrapper(gzip.open(names, "r")) as infile:
                f = infile.read()
                outfile.write(f)
            outfile.write("\n")
    return

if __name__ == "__main__":
    build_prokka_refgbff(sys.argv[1], sys.argv[2], sys.argv[3])