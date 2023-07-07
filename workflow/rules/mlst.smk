rule mlst:
    input:
        fna="data/interim/fasta/{strains}.fna",
    output:
        csv="data/interim/mlst/{strains}_ST.csv",
    conda:
        "../envs/mlst.yaml"
    log:
        "logs/mlst/mlst-{strains}.log",
    shell:
        """
        mlst --csv {input.fna} > {output.csv} 2>> {log}
        """
