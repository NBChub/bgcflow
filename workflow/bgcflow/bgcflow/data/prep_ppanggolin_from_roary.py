import sys
from pathlib import Path

import pandas as pd


def prep_ppanggolin_input(input_file, output_file):
    """
    Process a DataFrame from an single lined input separated by spaces,
    and output it as a tab separated list of file name and its corresponding path

    Parameters:
        input_file (str): Path to the input txt file.
        output_file (str): Path to the output file
    """
    # Read the input CSV file into a DataFrame
    df = pd.read_csv(input_file, sep=" ").T.reset_index(drop=False)

    # Perform the required transformations
    for i in df.index:
        data = Path(df.loc[i, "index"])
        genome_id = data.stem
        df.loc[i, "genome_id"] = genome_id

    df = df.set_index("genome_id")

    # Create the parent directory of the output file if it doesn't exist
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Export the modified DataFrame to the output file
    df.to_csv(output_file, sep="\t", header=False, columns=None)
    return


if __name__ == "__main__":
    prep_ppanggolin_input(sys.argv[1], sys.argv[2])
