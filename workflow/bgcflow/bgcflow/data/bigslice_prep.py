import pandas as pd
import sys


def prep_bigslice(df_gtdb_output, df_bigslice_output):
    df = pd.read_csv(df_gtdb_output)
    # Format to BiG-Slice
    df_bigslice = df.rename(columns={"genome_id": "#Genome Folder"})
    df_bigslice["#Genome Folder"] = df_bigslice["#Genome Folder"].apply(
        lambda x: f"{x}/"
    )
    df_bigslice.to_csv(df_bigslice_output, index=False, sep="\t")
    return None


if __name__ == "__main__":
    prep_bigslice(sys.argv[1], sys.argv[2])
