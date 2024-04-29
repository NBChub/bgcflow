import logging
import sys
from pathlib import Path

from Bio import SeqIO

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def extract_genbank_info(gb_file):
    """
    Extract information from a GenBank file and store it in a dictionary.

    Parameters:
    gb_file (str): The path to the GenBank file.

    Returns:
    dict: A dictionary containing the extracted information.
    """
    logging.info(f"Extracting information from {gb_file}")
    info_dict_all = {}
    for num, record in enumerate(SeqIO.parse(gb_file, "genbank")):
        source = [feat for feat in record.features if feat.type == "source"]
        assert (
            len(source) == 1
        ), f"Expected 1 source feature in the record, found {len(source)}"
        info_dict = {}
        info_dict["organism"] = source[0].qualifiers["strain"][0]
        info_dict["bases"] = len(record)
        info_dict["CDS"] = len([feat for feat in record.features if feat.type == "CDS"])
        info_dict["rRNA"] = len(
            [feat for feat in record.features if feat.type == "rRNA"]
        )
        info_dict["repeat_region"] = len(
            [feat for feat in record.features if feat.type == "repeat_region"]
        )
        info_dict["tRNA"] = len(
            [feat for feat in record.features if feat.type == "tRNA"]
        )
        info_dict["tmRNA"] = len(
            [feat for feat in record.features if feat.type == "tmRNA"]
        )
        info_dict_all[num] = info_dict
    return info_dict_all


def summarize_and_write(info_dict_all, output_file):
    """
    Summarize the information in a dictionary and write it to a text file.

    Parameters:
    info_dict_all (dict): The dictionary containing the information to summarize.
    output_file (str): The path to the output text file.

    Returns:
    None
    """
    logging.info(f"Summarizing information and writing to {output_file}")
    summary = {
        "organism": info_dict_all[0]["organism"],
        "contigs": len(info_dict_all),
        "bases": sum(info["bases"] for info in info_dict_all.values()),
        "CDS": sum(info["CDS"] for info in info_dict_all.values()),
        "rRNA": sum(info["rRNA"] for info in info_dict_all.values()),
        "repeat_region": sum(info["repeat_region"] for info in info_dict_all.values()),
        "tRNA": sum(info["tRNA"] for info in info_dict_all.values()),
        "tmRNA": sum(info["tmRNA"] for info in info_dict_all.values()),
    }

    with open(output_file, "w") as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")


def summarize_genbank(gb_file, output_file):
    """
    Extract information from a GenBank file, summarize it, and write the summary to a text file.

    Parameters:
    gb_file (str): The path to the GenBank file.
    output_file (str): The path to the output text file.

    Returns:
    None
    """
    logging.info(f"Starting to summarize GenBank file {gb_file}")
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    info_dict_all = extract_genbank_info(gb_file)
    logging.info(f"Writing summary to {output_file}")
    summarize_and_write(info_dict_all, output_file)


if __name__ == "__main__":
    summarize_genbank(sys.argv[1], sys.argv[2])
