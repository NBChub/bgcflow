rule mlst:
    input: 
        fna = "data/interim/fasta/{strains}.fna",
    output:
        csv = "data/interim/mlst/{strains}_ST.csv"
    conda:
        "../envs/mlst.yaml"
    shell:
        """
        mlst --csv {input.fna} > {output.csv}
        """