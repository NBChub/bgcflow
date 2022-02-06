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

    # Get taxonomic information
    df = pd.DataFrame([i for i in df_tax.gtdb_taxonomy])
    for col in df.columns:
        df[col] = df[col].apply(lambda x: x.split('__')[1])
    df["organism"] = df["species"]
    for idx in df.index:
        try:
            df.loc[idx, "species"] = df.loc[idx, "species"].split(' ')[1]
        except IndexError: # leave blank for empty taxonomy
            pass
    df.columns = [c.title() for c in df.columns]
    df.insert(0, "genome_id", df_tax.genome_id)

    # Get metadata
    for i in df_tax.index:
        try:
            metadata = df_tax.loc[i, "metadata"]
            for k in metadata.keys():
                if k.startswith("checkm"):
                    df.loc[i, k] = metadata[k]
                if k.startswith("ncbi"):
                    df.loc[i, k] = metadata[k]
        except:
            pass

    # save to file
    df.to_csv(df_gtdb_output, index=False)
    return None

if __name__ == "__main__":
    summarize_gtdb_json(sys.argv[1], sys.argv[2])