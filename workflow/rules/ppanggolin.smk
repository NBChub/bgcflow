
rule ppanggolin_prep:
    input:
        expand("data/interim/antismash/{strains}/", strains = STRAINS)
    output:
        table = "data/interim/antismash/ppanggolin_datasets.tsv",
    run:
        df = pd.DataFrame([[i,f"data/interim/antismash/{i}/{i}.gbk"] for i in df_samples.genome_id])
        df.to_csv(output.table, sep="\t", index=False, header=False)

rule ppanggolin:
    input:
        "data/interim/antismash/ppanggolin_datasets.tsv",
    output:
        directory("data/interim/ppanggolin/")
    conda:
        "../envs/ppanggolin.yaml"
    threads: 12
    shell:
        """
        ppanggolin workflow --anno {input} --cpu {threads} --output {output}
        """