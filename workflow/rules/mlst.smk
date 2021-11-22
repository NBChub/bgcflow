rule mlst:
    input: 
        gbk = "data/interim/prokka/{strains}/{strains}.gbk",
    output:
        csv = "data/interim/mlst/{strains}_ST.csv"
    conda:
        "../envs/mlst.yaml"
    shell:
        """
        mlst --help
        mlst --csv {input.gbk} > {output.csv}
        """