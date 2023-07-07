rule mash:
    input:
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        mash_infile="data/interim/mash/{name}/mash_in.txt",
        triangle_dist="data/interim/mash/{name}/triangle_distance_matrix.tsv",
    conda:
        "../envs/mash.yaml"
    threads: 32
    log:
        "logs/mash/mash-triangle-{name}.log",
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
        mash_matrix="data/interim/mash/{name}/triangle_distance_matrix.tsv",
    output:
        df_mash="data/processed/{name}/mash/df_mash.csv",
    log:
        "logs/mash/mash-convert-{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/convert_triangular_matrix.py {input.mash_matrix} {output.df_mash} 2>> {log}
        """
