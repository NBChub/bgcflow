import json
import logging
import sys
from pathlib import Path

from Bio import SeqIO

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def bgc_downstream_prep(input_dir, output_dir, selected_bgcs=False):
    """
    Given an antiSMASH directory, check for changed name
    """
    logging.info(f"Reading input directory: {input_dir}")
    path = Path(input_dir)
    if not path.is_dir():
        raise FileNotFoundError(f"No such file or directory: {path}")

    genome_id = path.name
    outpath = Path(output_dir) / genome_id
    outpath.mkdir(parents=True, exist_ok=True)
    logging.debug(f"Deducting genome id as {genome_id}")

    change_log = {genome_id: {}}
    ctr = 0
    matches = [Path(i).stem for i in selected_bgcs.split()]
    for gbk in path.glob("*.gbk"):
        if gbk.stem in matches:
            logging.debug(f"MATCH: {gbk.stem}")
            ctr = ctr + 1
            logging.info(f"Parsing file: {gbk.name}")
            region = SeqIO.parse(str(gbk), "genbank")
            for record in region:
                logging.info(f"{gbk} {record.id}")
                record_log = {}
                if "comment" in record.annotations:
                    filename = gbk.name
                    try:
                        original_id = record.annotations["structured_comment"][
                            "antiSMASH-Data"
                        ]["Original ID"].split()[0]
                    except KeyError:
                        original_id = record.id
                        logging.warning(
                            f"Found shortened record.id: {record.id} <- {original_id}."
                        )

                    # generate symlink
                    new_filename = filename.replace(record.id, original_id)
                    target_path = Path.cwd() / gbk  # target for symlink

                    link = outpath / new_filename

                    logging.info(f"Generating symlink: {link}")
                    try:
                        link.symlink_to(target_path)
                    except FileExistsError:
                        logging.warning(
                            f"Previous symlink exist, updating target: {link} -> {target_path}"
                        )
                        link.unlink()
                        link.symlink_to(target_path)

                    record_log["record_id"] = record.id
                    record_log["original_id"] = original_id
                    record_log["target_path"] = str(gbk)
                    record_log["symlink_path"] = str(link)
                else:
                    logging.warning(f"No Comments in record: {gbk.name}")

                change_log[genome_id][filename] = record_log
    # assert 1+1==3
    with open(
        outpath / f"{genome_id}-change_log.json", "w", encoding="utf8"
    ) as json_file:
        json.dump(change_log, json_file, indent=4)

    logging.info(f"{genome_id}: Job done!\n")
    return


if __name__ == "__main__":
    bgc_downstream_prep(sys.argv[1], sys.argv[2], sys.argv[3])
