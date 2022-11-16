import pandas as pd
from pathlib import Path
import json
import shutil
import logging
import sys

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def assess_gtdb_json_file(item):
    """
    Given a json file, assess whether the accession can be found via GTDB-API or not.
    Return a genome_id when accession cannot be found.
    """
    logging.info(f"Assessing {item}...")
    with open(item, "r") as json_file:
        data = json.load(json_file)
        genome_id = data["genome_id"]
        try:
            gtdb_release = data["gtdb_release"]
            metadata = data["metadata"]
            try:
                if type(metadata["genome"]["accession"]) == str:
                    logging.debug(
                        f"{genome_id} can be found via GTDB-API release {gtdb_release}"
                    )
                    return None
            except KeyError:
                if metadata["detail"] == "Genome not found":
                    logging.debug(
                        f"{genome_id} : {metadata['detail']} in GTDB-API release {gtdb_release}"
                    )
                    return genome_id
        except KeyError:
            logging.debug(f"{genome_id} does not have metadata")
            return genome_id


def generate_symlink_gtdbtk(input_fna, gtdb_json, outdir):
    """
    Given a gtdb_json file and an input_fna file, generate a symlink to a desired location
    if genome_id cannot be found via GTDB API
    """
    input_fna = Path(input_fna).resolve()
    gtdb_json = Path(gtdb_json).resolve()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    assert input_fna.is_file() and gtdb_json.is_file()

    genome_id = assess_gtdb_json_file(gtdb_json)
    if genome_id != None:
        outfile = outdir / f"{genome_id}.fna"
        logging.info(f"Generating input files for GTDB-tk: {outfile}")
        outfile.symlink_to(input_fna)
    return None


def gtdbtk_prep(fna_list, json_list, outdir):
    """
    Given a list of gtdb_json file and an list of fna, generate a symlinks to a desired location
    if genome_id cannot be found via GTDB API
    """
    shell_json_input = json_list.split()
    shell_fna_input = fna_list.split()
    for gtdb_json in shell_json_input:
        gid = Path(gtdb_json).stem
        input_fna = [fna for fna in shell_fna_input if gid in Path(fna).stem]
        generate_symlink_gtdbtk(str(input_fna[0]), str(gtdb_json), str(outdir))
    return


if __name__ == "__main__":
    gtdbtk_prep(sys.argv[1], sys.argv[2], sys.argv[3])
