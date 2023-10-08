import argparse

import pandas as pd


def convert_genome_to_species_mapping(input_file, output_file):
    """
    Convert the Genome_to_Species_Mapping.txt file from lsaBGC to the BGCFlow sample file format.

    Parameters:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    # output backbone
    df_samples = pd.DataFrame(
        columns=[
            "genome_id",
            "source",
            "organism",
            "genus",
            "species",
            "strain",
            "closest_placement_reference",
            "input_file",
        ]
    )

    # Read in the data
    df = pd.read_csv(input_file, sep="\t", header=None)

    for i in df.index:
        genus, species, rest = df.loc[i, 0].split("_", 2)
        genome_id = df.loc[i, 0].split("_", 3)[-1]
        genome_id = ".".join(genome_id.rsplit("_", 1))
        df_samples.loc[i, "genome_id"] = genome_id
        df_samples.loc[i, "source"] = "ncbi"
        df_samples.loc[i, "organism"] = " ".join([genus, species])
        df_samples.loc[i, "genus"] = genus
        df_samples.loc[i, "species"] = species
    df_samples.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Genome_to_Species_Mapping.txt file from lsaBGC to BGCFlow sample file format"
    )
    parser.add_argument("input_file", type=str, help="Path to the input file")
    parser.add_argument("output_file", type=str, help="Path to the output file")
    args = parser.parse_args()

    convert_genome_to_species_mapping(args.input_file, args.output_file)
