import logging
from pathlib import Path
import os
import sys

def create_symlink(source, destination):
    logging.info(f"Creating symlink from {source} to {destination}")
    source_path = Path(source).resolve()
    destination_path = Path(destination)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    if not destination_path.exists():
        os.symlink(source_path, destination)
        logging.info(f"Symlink created at {destination}")
    else:
        logging.warning(f"Destination {destination} already exists")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    source = sys.argv[1]
    destination = sys.argv[2]
    create_symlink(source, destination)