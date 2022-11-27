import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def bigscape_copy(input_index, output_main):
    input = Path(input_index)
    output = Path(output_main)
    output.mkdir(parents=True, exist_ok=True)
    for item in input.parent.glob("*"):
        target = output / item.name
        logging.debug(f"Generating symlink for: {target} --> {item.resolve()}")
        try:
            target.symlink_to(item.resolve(), target_is_directory=True)
        except FileExistsError as e:
            logging.debug(f"Got error:\n{e}")
    return


if __name__ == "__main__":
    bigscape_copy(sys.argv[1], sys.argv[2])
