import json
import os
import sys

import pandas as pd


def get_checkm_data(checkm_input_stats, checkm_json_folder, checkm_output_stats):
    """
    Read checkm output stats.tsv file, return json output for individual genomes and output stats table

    Parameters
    ----------
    1. checkm_input_stats : str / path
        Location of the output from checkm software 'storage/bin_stats_ext.tsv'
    2. checkm_json_folder : str / path
        Location of the json folder to store output for each genome
    3. checkm_output_stats : str / path
        Location of the output csv file in the processed/tables directory

    Returns
    -------
    1. json : .json file
        A json file for each genome summarizing the checkm output
    2. checkm_output_stats : str / path
        Output csv file in the processed/tables directory with checkm information for all genomes
    """

    # List of columns in df_checkm_stats.csv
    checkm_columns = [
        "Completeness",
        "Contamination",
        "GC",
        "GC std",
        "Genome size",
        "# ambiguous bases",
        "# scaffolds",
        "# contigs",
        "Longest scaffold",
        "Longest contig",
        "N50 (scaffolds)",
        "N50 (contigs)",
        "Mean scaffold length",
        "Mean contig length",
        "Coding density",
        "Translation table",
        "# predicted genes",
    ]

    df_checkm_input = pd.read_csv(
        checkm_input_stats, sep="\t", header=None, index_col=0
    )

    df_checkm_out = pd.DataFrame(index=df_checkm_input.index, columns=checkm_columns)
    df_checkm_out.index.name = "genome_id"

    for genome_id in df_checkm_out.index:
        json_genome_id = json.loads(df_checkm_input.loc[genome_id, 1].replace("'", '"'))
        json_path = os.path.join(checkm_json_folder, genome_id + ".json")

        with open(json_path, "w") as outfile:
            json.dump(json_genome_id, outfile)
        for col in df_checkm_out.columns:
            df_checkm_out.loc[genome_id, col] = json_genome_id[col]

    df_checkm_out.to_csv(checkm_output_stats)

    return df_checkm_out


if __name__ == "__main__":
    get_checkm_data(sys.argv[1], sys.argv[2], sys.argv[3])
