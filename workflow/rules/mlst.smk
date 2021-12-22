rule mlst:
    input: 
        fna = "data/interim/fasta/{strains}.fna",
    output:
        csv = "data/interim/mlst/{strains}_ST.csv"
    conda:
        "../envs/mlst.yaml"
    log: "workflow/report/logs/{strains}/mlst.log"
    shell:
        """
        mlst --csv {input.fna} > {output.csv}
        head {output.csv} >> {log}
        """