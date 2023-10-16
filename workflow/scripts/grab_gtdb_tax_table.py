import argparse
from pathlib import Path

import pandas as pd


def get_GTDB_tax_table(GTDB_Tax, outfile):
    """
        Download GTDB taxonomy data, process it, and save it as a TSV table to a specified file.

        Args:
            GTDB_Tax (str): URL to the GTDB taxonomy data file.
            outfile (str): Name of the output file to save the processed data.

        Returns:
            None

        Usage:
            python workflow/scripts/grab_gtdb_tax_table.py --url "https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv" --outfile config/Lactobacil
    lus_delbrueckii/bac120_taxonomy_r214.tsv
    """

    df = pd.read_csv(GTDB_Tax, sep="\t", header=None)

    # Add column names
    df.columns = ["user_genome", "classification"]

    # Modify the user_genome column
    df["user_genome"] = df["user_genome"].str.split("_", n=1).str[-1]

    # Print the updated DataFrame
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outfile, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Download and process GTDB taxonomy data."
    )
    parser.add_argument(
        "--url",
        required=True,
        help="URL to the GTDB taxonomy data file (e.g., https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv)",
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help="Name of the output TSV file to save the processed data",
    )
    args = parser.parse_args()

    get_GTDB_tax_table(args.url, args.outfile)


if __name__ == "__main__":
    main()
