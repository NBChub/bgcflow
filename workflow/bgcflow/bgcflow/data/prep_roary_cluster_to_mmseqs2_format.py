import sys
from pathlib import Path

import pandas as pd


def export_paired_values(input_file, output_file):
    """
    Read roary cluster output, create paired values, and export them to a text file.

    Parameters:
        input_file (str): Path to the input file.
        output_file (str): Path to the output text file.

    Usage example:
        input_file = "data/interim/roary/Lactobacillus_delbrueckii/clustered_proteins"
        output_file = 'paired_values.txt'
        export_paired_values(input_file, output_file)
    """
    # Read the input file into a DataFrame
    df_cluster = pd.read_csv(input_file, sep=":", header=None, index_col=0)

    # Process the DataFrame and create paired values
    data = {k: v.split("\t") for k, v in df_cluster.to_dict()[1].items()}
    paired_list = [(key, item) for key, values in data.items() for item in values]

    # Export the paired values to the output file
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as file:
        for pair in paired_list:
            file.write(f"{pair[0]}\t{pair[1]}\n")
    return


if __name__ == "__main__":
    export_paired_values(sys.argv[1], sys.argv[2])
