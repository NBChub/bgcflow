from pathlib import Path
from Bio import SeqIO
import logging
import json
import sys

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def bgc_downstream_prep(input_dir, output_dir): 
    """
    Given an antiSMASH directory, check for changed name
    """
    logging.info(f"Reading input directory: {input_dir}")
    path = Path(input_dir)

    if path.is_dir():
        project_type = "genome"
        pass
    elif path.is_file():
        assert path.suffix == ".gbk"
        path = path.parents[0]
        project_type = "bgc"
    else:
        raise FileNotFoundError(f'No such file or directory: {path}')

    genome_id = path.name
    logging.debug(f"Deducting genome id as {genome_id}")

    change_log = {genome_id : {}}
    original_change_log = {genome_id : {}}

    if project_type == "genome":
        iterate_over = '*.gbk'
    elif project_type == "bgc":
        iterate_over = Path(input_dir).name
    for gbk in path.glob(iterate_over):
        logging.info(f"Parsing file: {gbk.name}")
        region = SeqIO.parse(str(gbk), 'genbank')
        for record in region:
            record_log = {}
            original_record_log = {}
            if 'comment' in record.annotations:
                filename = gbk.name
                try:
                    original_id = record.annotations['structured_comment']['antiSMASH-Data']['Original ID'].split()[0]
                except KeyError:
                    original_id = record.id
                    logging.warning(f"Found shortened record.id: {record.id} <- {original_id}.")

                # generate symlink
                outpath = Path(output_dir) / genome_id
                outpath.mkdir(parents=True, exist_ok=True)

                new_filename = filename.replace(record.id, original_id)
                target_path = Path.cwd() / gbk # target for symlink

                link = outpath / new_filename

                logging.info(f'Generating symlink: {link}')
                try:
                    link.symlink_to(target_path)
                except FileExistsError as e:
                    logging.warning(f"Previous symlink exist, updating target: {link} -> {target_path}")
                    link.unlink()
                    link.symlink_to(target_path)

                record_log['record_id'] = record.id
                record_log['original_id'] = original_id
                record_log['target_path'] = str(gbk)
                record_log['symlink_path'] = str(link)

                original_record_log['record_id'] = record.id
                original_record_log['original_id'] = original_id
                original_record_log['target_path'] = str(gbk)

            change_log[genome_id][filename] = record_log
            original_change_log[genome_id][filename] = original_record_log

    with open(outpath / f"{genome_id}-change_log.json", 'w', encoding ='utf8') as json_file:
        json.dump(change_log, json_file, indent=4)

    # leave a change log in the original antismash parent folder
    logging.info("Adding new information to original record")
    original_log = Path(input_dir).parents[0] / f"{genome_id}.json"
    if original_log.is_file():
        pass
    else:
        with open(original_log, 'w', encoding ='utf8') as json_file:
            json.dump(original_change_log, json_file, indent=4)

    logging.info(f"{genome_id}: Job done!\n")
    return

if __name__ == "__main__":
    bgc_downstream_prep(sys.argv[1], sys.argv[2])