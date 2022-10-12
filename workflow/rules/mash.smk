rule mash:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        mash_infile = "data/interim/mash/{name}/mash_in.txt",
        triangle_dist = "data/interim/mash/{name}/triangle_distance_matrix.tsv"
    conda:
        "../envs/mash.yaml"
    threads: 32
    log: "workflow/report/logs/mash/mash-triangle-{name}.log"
    shell:
        """
        for fna in {input.fna} 
        do
            echo $fna >> {output.mash_infile}
        done
        (mash triangle -p {threads} -l {output.mash_infile} >> {output.triangle_dist}) 2>> {log}
        """

rule mash_convert:
    input:
        mash_matrix = "data/interim/mash/{name}/triangle_distance_matrix.tsv",
    output:
        df_mash = "data/processed/{name}/mash/df_mash.csv"
    run:
        df_raw = pd.read_csv(input.mash_matrix, index_col=0)
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
        df.to_csv(output.df_mash)
