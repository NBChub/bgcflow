import os
import sys
import pandas as pd

def fix_gtdb_taxonomy(samples_path, taxonomy_raw, taxonomy, meta):
    shell_input = samples_path.split()
    dfList = [pd.read_csv(s).set_index('genome_id', drop=False) for s in shell_input]
    df_samples = pd.concat(dfList, axis=0)

    df_tax = pd.read_csv(taxonomy_raw, sep="\t")
    df_tax = df_tax.set_index("#Genome folder", drop=False)

    dict_tax = {i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("custom")].genome_id}
    dict_tax.update({i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("patric")].genome_id})
    dict_tax.update({i : str(i)+"/" for i in df_samples[df_samples.source.eq("ncbi")].genome_id})

    container = []
    df_meta = pd.DataFrame(index=df_samples.index, columns=df_tax.columns)
    for i in df_samples.genome_id:
        line = df_tax.loc[dict_tax[i], :].to_dict()
        line.update({'#Genome folder': f"{i}/"})
        container.append(line)
        df_meta.loc[i, :] = df_tax.loc[dict_tax[i], :]
            
    df = pd.DataFrame(container)
    df.to_csv(taxonomy, sep="\t", index=False)
    df_meta.drop(columns='#Genome folder', inplace=True)
    df_meta.to_csv(meta)
    return None

if __name__ == "__main__":
    fix_gtdb_taxonomy(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])