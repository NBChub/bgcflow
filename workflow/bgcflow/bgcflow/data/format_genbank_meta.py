import json
import logging
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def get_git_version():
    """
    Get the sha1 of the current git version
    """

    git_version = ""
    try:
        version_cmd = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        status_cmd = subprocess.run(
            ["git", "status", "--porcelain"],
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        git_version = str(version_cmd.stdout.strip())
        changes = str(status_cmd.stdout).strip()
        if changes:
            git_version += "(changed)"
    except OSError:
        pass

    return git_version


def get_version(version):
    """
    Get the current version string
    """

    git_version = get_git_version()
    if git_version:
        version += "-%s" % git_version
    return version


def add_bgcflow_comments(gbk_in_path, version, json_path, genome_id, gbk_out_path):
    """
    Add bgcflow meta-annotation to genbank output
    """

    version = get_version(version)

    logging.info(f"Formatting genbank file: {gbk_in_path}")

    temp_gbk = Path(gbk_in_path).parent / f"{genome_id}-change_log.gbk"
    try:
        [record for record in SeqIO.parse(gbk_in_path, "genbank")]
    except ValueError as e:
        logging.warning(f"Parsing fail: {e}")
        logging.info("Attempting to fix genbank file...")
        change_log_path = Path(gbk_in_path).parent / f"{genome_id}-change_log.json"
        change_log = correct_collided_headers(
            gbk_in_path, genome_id, temp_gbk, json_dump=change_log_path
        )
        logging.info("Retry parsing with Bio.SeqIO...")

    if temp_gbk.is_file():
        records = SeqIO.parse(temp_gbk, "genbank")
        temp_gbk.unlink()
    else:
        records = SeqIO.parse(gbk_in_path, "genbank")

    df_gtdb_meta = pd.read_json(json_path, orient="index").T
    df_gtdb_meta.fillna("Unclassified", inplace=True)
    df_tax = pd.DataFrame([i for i in df_gtdb_meta.gtdb_taxonomy])
    for col in df_tax.columns:
        df_tax[col] = df_tax[col].apply(lambda x: x.split("__")[1])
    df_tax.columns = [c.title() for c in df_tax.columns]
    tax_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    taxonomy_str = df_tax.loc[0, tax_levels].tolist()
    logging.debug(f"Taxonomy found: {taxonomy_str}")

    bgcflow_comment = (
        "##BGCflow-Data-START##\n"
        "Version      :: {version}\n"
        "Run date     :: {date}\n"
    ).format(
        version=version,
        date=str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
    )

    new_records = []
    ctr = 0
    for record in records:
        logging.info(f"Formatting record {ctr}...")
        comment = bgcflow_comment
        if "comment" in record.annotations:
            record.annotations["comment"] += "\n" + comment
        else:
            record.annotations["comment"] = comment

        try:
            if not change_log[ctr]["accept_format"]:
                record.annotations[
                    "comment"
                ] += f"Original ID  :: {change_log[ctr]['original_id']}\n"
        except UnboundLocalError:
            pass

        record.annotations["comment"] += "##BGCflow-Data-END##"

        if "organism" in record.annotations:
            organism = record.annotations["organism"]
            if "Unclassified" in organism:
                record.annotations["organism"] = organism.split(" Unclassified")[
                    0
                ].strip()

        record.annotations["taxonomy"] = taxonomy_str

        new_records.append(record)

        ctr = ctr + 1

    logging.info(f"Writing final output: {gbk_out_path}")
    with open(gbk_out_path, "w") as output_handle:
        SeqIO.write(new_records, output_handle, "genbank")


# -- Correcting IDs --#
# 1. Traverse text, make a library of record header, length, and ids
# 2. Check the library for collided locus name and length
# 3. Propose changes
# 4. Make sure changes is unique
# 5. Write a new genbank file and recorded changes as json file


def correct_collided_headers(
    gbk_in_path, accession_id, outfile, json_dump="header.json"
):
    record_headers = get_record_headers(gbk_in_path)
    record_headers = shorten_record_headers(record_headers, accession_id)

    with open(json_dump, "w", encoding="utf8") as json_file:
        json.dump(record_headers, json_file, indent=4)

    logging.info(f"Writing result to: {outfile}")
    with open(gbk_in_path) as f:
        s = f.read()

    with open(outfile, "w") as f:
        for k in record_headers.keys():
            if not record_headers[k]["accept_format"]:
                s = s.replace(
                    record_headers[k]["original_header"],
                    record_headers[k]["new_header"],
                )
        f.write(s)
    return record_headers


def get_record_headers(gbk_in_path):
    """
    Get the header information for each records in a genbank file.
    Returns a json-like python dictionary.
    """
    record_headers = {}

    # Find all headers
    with open(gbk_in_path, "r") as file:
        logging.info(f"Reading file as text: {gbk_in_path}")
        ctr = 0
        for line in file:
            if line.startswith("LOCUS"):
                record_headers[ctr] = {"original_header": line}
            if "source" in line:
                length = line.split()[-1].split("..")[-1]
                record_headers[ctr]["length"] = length
                ctr = ctr + 1
    logging.debug(f"Found {len(record_headers)} records.")

    # Check for collided headers
    logging.info("Checking header format...")
    for k in record_headers.keys():
        query = record_headers[k]["original_header"].split()
        if query != 7 and query[3] != "bp":
            logging.warning("Record name and length collide in the LOCUS line")
            query_length = record_headers[k]["length"]
            if query_length in query[1]:
                original_id = query[1].replace(query_length, "")
                record_headers[k]["original_id"] = original_id
                record_headers[k]["accept_format"] = False
        else:
            record_headers[k]["original_id"] = query[1]
            record_headers[k]["accept_format"] = True
    return record_headers


def modify_id(accession_id, original_id, locus_index, record_counts, id_collection):
    """
    Attempt to fix collided name and length in the locus name by shortening id
    """
    logging.info(f"{original_id}: Attempting to change locus name...")

    # For refseq record, remove the accession prefix (first two digits) from string
    if accession_id.startswith("GCF"):
        logging.debug(
            f"{original_id}: Assuming locus came from Refseq. Removing refseq accession prefix."
        )
        refseq_type, refseq_number, refseq_version = re.split("_|[.]", original_id)
        new_id = f"{refseq_number}.{refseq_version}"

    elif accession_id.startswith("GCA"):
        logging.info(
            f"{original_id}: Assuming locus came from Genbank. Removing version from locus name..."
        )
        new_id, genbank_version = re.split("[.]", original_id)

    # For unknown source remove last 4 digit, add the contig index in front.
    # example: NZAJABAQG010000001.1 --> c1|NZAJABAQG010000
    else:
        logging.info(
            f"{original_id}: Cannot determine source. Shortening locus name..."
        )
        digit = len(str(record_counts))
        contig_number = str(locus_index + 1)

        new_id = f"c{contig_number.zfill(digit)}_{original_id[:-(digit + 4)]}"

    # make sure new_id is unique
    logging.info(f"{original_id}: Making sure id is unique...")
    if new_id in id_collection:
        raise
    else:
        logging.debug(f"{original_id}: shortened to {new_id}")
    return new_id


def shorten_record_headers(record_headers, accession_id):
    """
    Shorten record headers
    """
    record_counts = len(record_headers.keys())
    id_collection = [record_headers[k]["original_id"] for k in record_headers.keys()]

    for k in record_headers.keys():
        if record_headers[k]["accept_format"]:
            pass
        else:
            original_id = record_headers[k]["original_id"]
            old_header = record_headers[k]["original_header"]

            # Correct header name
            new_id = modify_id(
                accession_id, original_id, k, record_counts, id_collection
            )
            new_id_replacement_value = (
                f"{new_id}{(len(original_id) - len(new_id)) * ' '}"
            )
            new_header = old_header.replace(original_id, new_id_replacement_value)

            # Add value to dict
            record_headers[k]["new_id"] = new_id
            record_headers[k]["new_header"] = new_header

    return record_headers


if __name__ == "__main__":
    add_bgcflow_comments(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    )
