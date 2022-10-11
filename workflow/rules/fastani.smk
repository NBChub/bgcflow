rule fastani:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fastani_infile = "data/interim/fastani/{name}/fastani_in.txt",
        fastani_out = "data/interim/fastani/{name}/fastani_out.tsv",
        fastani_matrix = "data/interim/fastani/{name}/fastani_out.tsv.matrix",
    conda:
        "../envs/fastani.yaml"
    threads: 32
    log: "workflow/report/logs/fastani/fastani-{name}.log"
    shell:
        """
        for fna in {input.fna} 
        do
            echo $fna >> {output.fastani_infile}
        done
        fastANI --ql {output.fastani_infile} --rl {output.fastani_infile} -t {threads} --matrix -o {output.fastani_out} 2>> {log} 
        """

rule fastani_convert:
    input:
        fastani_matrix = "data/interim/fastani/{name}/fastani_out.tsv.matrix",
    output:
        df_fastani = "data/processed/{name}/fastani/df_fastani.csv"
    run:
        df_raw = pd.read_csv(input.fastani_matrix, index_col=0)
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
        df.to_csv(output.df_fastani)