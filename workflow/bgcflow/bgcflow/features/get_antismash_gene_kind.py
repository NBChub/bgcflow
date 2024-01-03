import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)


def extract_gene_kinds(input_path, output_file):
    """
    Extracts the 'gene_kind' qualifier from each CDS feature in a GenBank file or all GenBank files in a directory,
    and writes the results to a CSV file.

    Parameters:
    input_path (str): The path to the GenBank file or directory containing GenBank files.
    output_file (str): The path to the output CSV file.
    """

    logging.info(f"Processing input: {input_path}")

    output = {}

    # Create a Path object from the input path
    input_path = Path(input_path)

    # Check if the input path is a directory
    if input_path.is_dir():
        # Iterate over all GenBank files in the directory
        for genbank_file in input_path.glob("*.gbk"):
            logging.info(f"Processing GenBank file: {genbank_file}")
            process_genbank_file(genbank_file, output)
    else:
        # Process the single GenBank file
        logging.info(f"Processing GenBank file: {input_path}")
        process_genbank_file(input_path, output)

    # Create the parent directory of the output file if it doesn't exist
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Write the results to the output CSV file
    pd.DataFrame.from_dict({"gene_kind": output}).to_csv(output_file, header=False)

    logging.info(f"Results written to: {output_file}")


def process_genbank_file(genbank_file, output):
    """
    Processes a GenBank file and extracts the 'gene_kind' qualifier from each CDS feature.

    Parameters:
    genbank_file (Path): The path to the GenBank file.
    output (dict): The dictionary to store the results.
    """

    # Open the GenBank file
    with open(genbank_file, "r") as input_handle:
        # Parse the GenBank file
        for record in SeqIO.parse(input_handle, "genbank"):
            # Iterate over each feature in the record
            for feature in record.features:
                # Check if the feature is a CDS and has the 'gene_kind' qualifier
                if feature.type == "CDS" and "gene_kind" in feature.qualifiers:
                    if "gene_kind" not in feature.qualifiers:
                        gene_kind = feature.qualifiers["gene_kind"][0]
                    else:
                        logging.warning(
                            f"Feature {feature} does not have a 'gene_kind' qualifier"
                        )
                        gene_kind = "other"
                    if "locus_tag" in feature.qualifiers:
                        locus_tag = feature.qualifiers["locus_tag"][0]
                    elif "gene" in feature.qualifiers:
                        locus_tag = feature.qualifiers["gene"][0]
                    else:
                        locus_tag = feature.qualifiers["protein_id"][0]
                    output[locus_tag] = gene_kind


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Extract gene kinds from a GenBank file or directory"
    )

    # Add the arguments
    parser.add_argument(
        "input_path",
        metavar="input_path",
        type=str,
        help="the path to the GenBank file or directory",
    )
    parser.add_argument(
        "output_file",
        metavar="output_file",
        type=str,
        help="the path to the output CSV file",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the arguments
    extract_gene_kinds(args.input_path, args.output_file)
