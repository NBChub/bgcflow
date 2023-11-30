import logging
import sys

from BCBio import GFF
from Bio import SeqIO

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def convert_gbk_to_gff3(input_file, output_file, write_fasta=True):
    """
    Convert a GenBank (.gbk) file to a GFF3 file with sequence.

    Parameters:
    input_file (str): Path to the input GenBank file.
    output_file (str): Path to the output GFF3 file.
    """
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        genbank = SeqIO.parse(input_handle, "genbank")
        GFF.write(genbank, output_handle, include_fasta=write_fasta)


if __name__ == "__main__":
    write_fasta = sys.argv[3].lower() == "true"
    convert_gbk_to_gff3(sys.argv[1], sys.argv[2], write_fasta)
