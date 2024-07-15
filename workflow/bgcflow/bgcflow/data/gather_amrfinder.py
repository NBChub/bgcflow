import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def cleanup_amrfinder(input_file, genome_id=None):
    """
    Cleans up the AMRFinder output file by filling null values in the 'Protein identifier' column
    with a combination of 'Contig id', 'Start', and 'Stop' coordinates. Adds a 'genome_id' column.

    Parameters:
    input_file (str): Path to the input file.
    genome_id (str, optional): Genome identifier. If not provided, it will be derived from the input file name.

    Returns:
    pd.DataFrame: Cleaned DataFrame.
    """
    input_path = Path(input_file)
    if genome_id is None:
        genome_id = input_path.stem
    df = pd.read_csv(input_path, sep="\t")

    # Fill null values in 'Protein identifier' with a combination of 'Contig id', 'Start', and 'Stop'
    df["Protein identifier"] = df.apply(
        lambda row: row["Protein identifier"]
        if pd.notnull(row["Protein identifier"])
        else f"{row['Contig id']}_{row['Start']}_{row['Stop']}",
        axis=1,
    )
    df["genome_id"] = genome_id

    return df


def gather_amrfinder(input_list, output_file):
    """
    Gathers and cleans up multiple AMRFinder output files listed in an input file,
    concatenates them into a single DataFrame, and writes the result to an output file.

    Parameters:
    input_list (str): Path to the file containing a list of input file paths.
    output_file (str): Path to the output file where the concatenated DataFrame will be saved.

    Returns:
    None
    """
    dataframes = []

    with open(input_list, "r") as f:
        data = f.read()

    dataframes = pd.concat([cleanup_amrfinder(i) for i in data.split()])
    dataframes.to_csv(output_file, index=False)
    return


if __name__ == "__main__":
    gather_amrfinder(sys.argv[1], sys.argv[2])
