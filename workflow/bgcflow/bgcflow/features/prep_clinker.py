import logging
import sys
from pathlib import Path

from Bio import SeqIO


def correct_gbks_for_mmseqs2(filepath):
    with open(filepath, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            try:
                correct_id = record.annotations["structured_comment"]["antiSMASH-Data"][
                    "Original ID"
                ].split()[0]
            except KeyError:
                logging.info(record.annotations["structured_comment"]["antiSMASH-Data"])
                correct_id = record.id
            region = filepath.stem.split(".")[-1]
            record.id = f"{correct_id}.{region}"
            record.name = f"{correct_id}.{region}"
    return record, record.name


def prepare_gbks_for_mmseqs2(input_files, output_folder):
    """
    Correct shortened names and merge genbank together
    """
    # Define the input files and output file name
    input_files = input_files.split()

    # Create a list to hold the SeqRecord objects
    records = {}

    # Loop through each input file and append the records to the list
    for file in input_files:
        data, name = correct_gbks_for_mmseqs2(Path(file))
        records[name] = data

    # Write out the concatenated records as a new GenBank file
    for k, v in records.items():
        outfile = Path(output_folder) / f"{k}.gbk"
        outfile.parent.mkdir(exist_ok=True, parents=True)
        with open(outfile, "w") as f:
            SeqIO.write(v, f, "genbank")

    return


if __name__ == "__main__":
    prepare_gbks_for_mmseqs2(sys.argv[1], sys.argv[2])
