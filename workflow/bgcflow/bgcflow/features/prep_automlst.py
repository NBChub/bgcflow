import sys

from Bio import SeqIO


def update_gbk_automlst(input_gbk, out_gbk, genome_id):
    """
    Update organism field with accession ID for autoMLST run
    """

    input_handle = open(input_gbk, "r")
    new_seq_records = []

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        seq_record.annotations["organism"] = genome_id
        new_seq_records.append(seq_record)

    SeqIO.write(new_seq_records, out_gbk, "genbank")
    input_handle.close()

    return None


if __name__ == "__main__":
    update_gbk_automlst(sys.argv[1], sys.argv[2], sys.argv[3])
