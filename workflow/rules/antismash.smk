rule antismash_db_setup:
    output:
        directory("resources/antismash_db"),
    conda:
        "../envs/antismash.yaml"
    shell:  
        """
        download-antismash-databases --database-dir resources/antismash_db
        """

rule antismash:
    input: 
        gbk = "data/processed/genbank/{strains}.gbff",
        resources = "resources/antismash_db/"
    output:
        folder = directory("data/interim/antismash/{strains}"),
        gbk = "data/interim/antismash/{strains}/{strains}.gbk"
    conda:
        "../envs/antismash.yaml"
    threads: 8
    shell:
        """
        antismash --genefinding-tool prodigal --output-dir {output.folder} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} -v
        """