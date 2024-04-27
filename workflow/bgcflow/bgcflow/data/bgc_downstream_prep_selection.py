import json
import logging
import sys
from pathlib import Path

from Bio import SeqIO

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def generate_symlink(path, genome_id, output_dir, selected_bgcs=False):
    """
    Given an antiSMASH directory, check for changed name
    """
    outpath = Path(output_dir) / genome_id
    outpath.mkdir(parents=True, exist_ok=True)
    logging.debug(f"Deducting genome id as {genome_id}")
    ctr = 0
    matches = selected_bgcs.stem
    for gbk in path.glob("*.gbk"):
        if gbk.stem in matches:
            logging.debug(f"Found match: {gbk.stem}")
            filename = gbk.name
            ctr = ctr + 1
            logging.info(f"Parsing file: {gbk.name}")
            region = SeqIO.parse(str(gbk), "genbank")
            for record in region:
                logging.debug(f"Processing: {gbk.name}: {record.id}")
                record_log = {}
                if "structured_comment" in record.annotations:
                    try:
                        original_id = record.annotations["structured_comment"][
                            "antiSMASH-Data"
                        ]["Original ID"].split()[0]
                    except KeyError:
                        original_id = record.id
                        logging.warning(
                            f"Found shortened record.id: {record.id} <- {original_id}."
                        )
                else:
                    raise ValueError(f"No Structured Comments in record: {gbk.name}")

                if (":" in str(record.description)) or (":" in original_id):
                    logging.warning(
                        f"Illegal character ':' found in genbank description, removing: {record.description}"
                    )
                    # Remove colon from description
                    record.description = record.description.replace(":", "")
                    original_id = original_id.replace(":", "")

                    # Rename antiSMASH comment
                    if "structured_comment" in record.annotations:
                        if (
                            "Original ID"
                            in record.annotations["structured_comment"][
                                "antiSMASH-Data"
                            ]
                        ):
                            record.annotations["structured_comment"]["antiSMASH-Data"][
                                "Original ID"
                            ] = original_id

                    # Write new GenBank file
                    new_filename = filename.replace(record.id, original_id)
                    with open(outpath / new_filename, "w") as output_handle:
                        SeqIO.write(record, output_handle, "genbank")
                    link = outpath / new_filename
                else:
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

                    # Assert that the symlink was correctly generated
                    assert link.is_symlink(), f"Failed to create symlink: {link}"
                    assert (
                        link.resolve() == target_path.resolve()
                    ), f"Symlink {link} does not point to the correct target: {target_path}"

                record_log["record_id"] = record.id
                record_log["original_id"] = original_id
                record_log["target_path"] = str(gbk)
                record_log["symlink_path"] = str(link)

                change_log = {filename: record_log}
    return change_log


def bgc_downstream_prep(input_file, output_dir):
    logging.info(f"Reading input file: {input_file}")
    with open(input_file, "r") as file:
        file_paths = [Path(f) for f in file.read().splitlines()]
    change_log_containers = {}
    for num, selected_bgcs in enumerate(file_paths):
        input_dir = selected_bgcs.parent
        logging.info(f"Reading input directory: {input_dir}")
        path = Path(input_dir)
        if not path.is_dir():
            raise FileNotFoundError(f"No such file or directory: {path}")

        # check if it has complete antiSMASH results
        if (path / f"{path.name}.json").is_file():
            logging.info("Found full antiSMASH record")
            genome_id = path.name
        else:
            logging.warning("No full antiSMASH record found, unknown genome id")
            genome_id = "unknown_genome_id"

        assert selected_bgcs.exists(), f"File does not exist: {selected_bgcs}"
        region_change_log = generate_symlink(path, genome_id, output_dir, selected_bgcs)
        change_log_containers[num] = {
            "genome_id": genome_id,
            "value": region_change_log,
        }
    change_logs = {}
    genome_ids = set(v["genome_id"] for v in change_log_containers.values())
    for genome_id in genome_ids:
        change_log = {}
        for v in change_log_containers.values():
            if v["genome_id"] == genome_id:
                entry_name = list(v["value"].keys())[0]
                change_log[entry_name] = v["value"][entry_name]
        change_logs[genome_id] = change_log
    logging.debug(change_logs)

    for genome_id in change_logs.keys():
        outpath = Path(output_dir) / genome_id
        with open(
            outpath / f"{genome_id}-change_log.json", "w", encoding="utf8"
        ) as json_file:
            json.dump({genome_id: change_logs[genome_id]}, json_file, indent=4)
        logging.info(f"{genome_id}: Job done!\n")


if __name__ == "__main__":
    bgc_downstream_prep(sys.argv[1], sys.argv[2])
