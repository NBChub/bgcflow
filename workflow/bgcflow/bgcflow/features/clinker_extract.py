import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def extract_clinker_to_table(clinker_txt, cutoff=0.7):
    """
    Convert clinker output to table format

    Arguments:
        - clinker_txt -- clinker output file
        - cutoff -- cutoff for identity (default: 0.7)

    Returns:
        - df_clinker -- dataframe of clinker output
    """
    with open(clinker_txt, "r") as f:
        clinker_data = f.readlines()

    clinker_dict = {}
    container = {}
    for num, line in enumerate(clinker_data):
        if line.startswith("----"):
            header = clinker_data[num - 1].strip("\n")
            if header not in clinker_dict.keys():
                clinker_dict[header] = []
        else:
            if line.startswith("----") or line.startswith("Query") or "vs" in line:
                pass
            else:
                data = line.strip("\n").split()
                seq_id = header.split(" vs ")
                if len(data) == 0:
                    data = [None, None, None, None]
                data = seq_id + data
                assert len(data) == 6, data
                clinker_dict[header].append(data)
                container[num] = data

    df_clinker = pd.DataFrame.from_dict(
        container,
        orient="index",
        columns=[
            "seq_id",
            "seq_id2",
            "locus_tag",
            "locus_tag2",
            "identity",
            "similarity",
        ],
    )
    df_clinker.identity = df_clinker.identity.astype(float)
    df_clinker = df_clinker[df_clinker.identity >= cutoff]
    return df_clinker


def build_cds_dict(gbk_dir):
    """ """
    gbk_dir = Path(gbk_dir)
    gbk_files = gbk_dir.glob("*.gbk")

    cds_dict = {}
    for gbk in gbk_files:
        seq_id = gbk.name
        logging.info(f"Processing {seq_id}")
        with open(gbk, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        translation = feature.qualifiers["translation"][0]
                        strand = feature.location.strand
                        if strand > 0:
                            start, end = int(feature.location.start), int(
                                feature.location.end
                            )
                        elif strand < 0:
                            end, start = int(feature.location.start), int(
                                feature.location.end
                            )
                        try:
                            locus_tag = feature.qualifiers["locus_tag"][0]
                        except KeyError:
                            logging.warning(f"{seq_id}: {feature.qualifiers.keys()}")
                            locus_tag = feature.qualifiers["protein_id"][0]
                            logging.warning(locus_tag)
                        if seq_id.startswith("BGC"):
                            locus_tag = feature.qualifiers["protein_id"][0]
                        cds_dict[locus_tag] = {
                            "seq_id": seq_id,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "translation": translation,
                        }
    return cds_dict


def annotate_clinker_table(df_clinker, cds_dict, outfile):
    """
    Annotate clinker table with coding sequence information

    Arguments:
        - df_clinker -- dataframe of clinker output
        - cds_dict -- dictionary of coding sequence information
        - outfile -- output file
    """
    for i in df_clinker.index:
        locus_tag1 = df_clinker.loc[i, "locus_tag"]
        locus_tag2 = df_clinker.loc[i, "locus_tag2"]
        try:
            if locus_tag1 is not None:
                df_clinker.loc[i, "start"] = int(cds_dict[locus_tag1]["start"])
                df_clinker.loc[i, "end"] = int(cds_dict[locus_tag1]["end"])
        except KeyError as e:
            logging.warning(f"{e}, {df_clinker.loc[i, 'seq_id']}")
        try:
            if locus_tag2 is not None:
                df_clinker.loc[i, "start2"] = int(cds_dict[locus_tag2]["start"])
                df_clinker.loc[i, "end2"] = int(cds_dict[locus_tag2]["end"])
        except KeyError as e:
            logging.warning(f"{e}, {df_clinker.loc[i, 'seq_id2']}")

    df_clinker.to_csv(outfile, index=False)
    return


def clinker_extract(clinker_txt, gbk_dir, outfile, cutoff=0.7):
    cutoff = float(cutoff)
    df_clinker = extract_clinker_to_table(clinker_txt, cutoff=cutoff)
    cds_dict = build_cds_dict(gbk_dir)
    annotate_clinker_table(df_clinker, cds_dict, outfile)
    return


if __name__ == "__main__":
    clinker_extract(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
