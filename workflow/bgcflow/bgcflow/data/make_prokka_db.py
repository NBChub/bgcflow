import ncbi_genome_download as ngd
import pandas as pd
import os, sys, glob, gzip, io
from pathlib import Path
import logging

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.INFO)

def build_prokka_refgbff(name, prokka_db_table, outfile):
    """
    Given a project name and csv table path, download all NCBI accession id
    and returns a concatenated gbk file for prokka --proteins params.
    """
    download_dir = Path(f"resources/prokka_db/{name}")
    prokka_db_table = Path(prokka_db_table).resolve()

    # Generate download directory
    download_dir.mkdir(parents=True, exist_ok=True)

    # Get lists of accession ids and its sources
    logging.info(f"Reading input file: {prokka_db_table}...")
    if prokka_db_table.suffix == ".csv":
        df = pd.read_csv(prokka_db_table)
    elif prokka_db_table.suffix == ".json":
        df = pd.read_json(prokka_db_table)

    ngd_input = {"refseq" : [],
                 "genbank" : []
                }

    logging.info("Donwloading genomes from NCBI...")
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
            logging.debug(f"Downloading {s} genome: {acc}")
            ngd.download(section=s,
                         file_formats='genbank', 
                         assembly_accessions=acc_list,  
                         output=download_dir, 
                         groups="bacteria")
    
    # Concatenate gbff
    logging.info("Concatenating genomes into a single file...")
    reference_files = glob.glob(f"{download_dir}/*/*/*/*.gbff.gz")

    with open(outfile, 'w') as f_out:
        for names in reference_files:
            with io.TextIOWrapper(gzip.open(names, "r")) as infile:
                f = infile.read()
                f_out.write(f)
            f_out.write("\n")

        logging.info(f"Reference genomes saved as: {str(f_out)}")

    return

if __name__ == "__main__":
    build_prokka_refgbff(sys.argv[1], sys.argv[2], sys.argv[3])