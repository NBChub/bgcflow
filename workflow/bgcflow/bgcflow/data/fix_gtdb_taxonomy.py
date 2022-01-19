import os, sys, json
import pandas as pd

def summarize_gtdb_json(accession_list, df_gtdb_output):
    # read json
    accession = accession_list.split()
    out = []
    for a in accession:
        with open(a, "r") as f:
            out.append(json.load(f))
    df_tax = pd.DataFrame(out)

    # Formatting
    df_tax = pd.concat([df_tax, pd.DataFrame([i for i in df_tax.gtdb_taxonomy])], axis=1)
    df = pd.DataFrame([i for i in df_tax.gtdb_taxonomy])
    for col in df.columns:
        df[col] = df[col].apply(lambda x: x.split('__')[1])
    df["organism"] = df["species"]
    df["species"] = df["species"].apply(lambda x: x.split(' ')[1])
    df.columns = [c.title() for c in df.columns]
    df.insert(0, "genome_id", df_tax.genome_id)
    df.to_csv(df_gtdb_output, index=False)
    return None

if __name__ == "__main__":
    summarize_gtdb_json(sys.argv[1], sys.argv[2])