import json
import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def assess_gtdb_json_file(item):
    """
    Given a json file, assess whether the accession can be found via GTDB-API or not.
    Return a genome_id when accession cannot be found.
    """
    logging.info(f"Assessing {item}")
    with open(item, "r") as json_file:
        data = json.load(json_file)
        genome_id = data["genome_id"]
        try:
            gtdb_release = data["gtdb_release"]
            metadata = data["metadata"]
            if "Genome not found" in metadata["detail"]:
                logging.debug(
                    f" - {genome_id} : {metadata['detail']} in GTDB-API release {gtdb_release}"
                )
                return genome_id
            elif type(metadata["genome"]["accession"]) == str:
                logging.debug(
                    f" - {genome_id} can be found via GTDB-API release {gtdb_release}"
                )
                return None

        except KeyError:
            logging.debug(f" - {genome_id} does not have metadata")
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
    if genome_id is not None:
        outfile = outdir / f"{genome_id}.fna"
        logging.info(f"Generating input files for GTDB-tk: {outfile}")
        outfile.symlink_to(input_fna)
        return None
    else:
        return genome_id


def input_handling(input_list, category, suffix=".json"):
    input_list = Path(input_list)
    if input_list.is_file() and input_list.suffix == suffix:
        logging.info(f"Getting {category} from a single file: {input_list}")
        input_list_files = input_list

    elif input_list.is_file() and input_list.suffix == ".txt":
        logging.info(f"Getting {category} from a text file: {input_list}")
        with open(input_list, "r") as file:
            file_content = [i.strip("\n") for i in file.readlines()]
            if len(file_content) == 1:
                # Paths space-separated on a single line
                logging.info(
                    " - Detecting space-separated input in a single line format."
                )
                paths = file_content[0].split()
            else:
                # Paths written on separate lines
                logging.info(" - Detecting input in a multi-line format.")
                paths = file_content
            input_list_files = [
                Path(path) for path in paths if Path(path).suffix == suffix
            ]
    else:
        input_list_files = [
            Path(file)
            for file in str(input_list).split()
            if Path(file).suffix == suffix
        ]
        logging.info(
            f"Getting {category} from the given list of {len(input_list_files)} files..."
        )
    return input_list_files


def gtdbtk_prep(fna_list, json_list, outdir, output_txt):
    """
    Given a list of gtdb_json file and an list of fna, generate a symlinks to a desired location
    if genome_id cannot be found via GTDB API
    """
    shell_json_input = input_handling(json_list, "taxonomy json")
    shell_fna_input = input_handling(fna_list, "fna files", suffix=".fna")

    input_list = []
    for gtdb_json in shell_json_input:
        gid = Path(gtdb_json).stem
        input_fna = [fna for fna in shell_fna_input if gid in Path(fna).stem]
        logging.info(
            f"Found {gid} in {[str(i) for i in input_fna]}. Generating symlink..."
        )
        fnafile = generate_symlink_gtdbtk(
            str(input_fna[0]), str(gtdb_json), str(outdir)
        )
        if fnafile is not None:
            input_list.append(fnafile)
    with open(output_txt, "w") as f:
        for item in input_list:
            f.write("%s\n" % item)
    return


if __name__ == "__main__":
    gtdbtk_prep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
