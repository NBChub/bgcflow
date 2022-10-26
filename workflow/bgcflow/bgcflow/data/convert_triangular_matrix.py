import pandas as pd
import sys

def triangle_convert(matrix, outfile):
    df_raw = pd.read_csv(matrix, index_col=0)
    genome_id_list = []
    for idx in df_raw.index:
        genome_id = idx.split('\t')[0].split('/')[-1].split('.fna')[0]
        genome_id_list.append(genome_id)

    df = pd.DataFrame(0, index=genome_id_list, columns=genome_id_list)
    for idx in df_raw.index:
        genome_id = idx.split('\t')[0].split('/')[-1].split('.fna')[0]
        for cntr in range(len(idx.split('\t'))):
            if cntr > 0:
                df.loc[genome_id, genome_id_list[cntr-1]] = float(idx.split('\t')[cntr])
                df.loc[genome_id_list[cntr-1], genome_id] = float(idx.split('\t')[cntr])
    df.index.name = 'genome_id'
    df.to_csv(outfile)
    return

if __name__ == "__main__":
    triangle_convert(sys.argv[1], sys.argv[2])