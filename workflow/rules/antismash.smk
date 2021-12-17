rule antismash_db_setup:
    output:
        directory("resources/antismash_db"),
    conda:
        "../envs/antismash.yaml"
    log: "workflow/report/antismash_db_setup.log"
    shell:  
        """
        download-antismash-databases --database-dir resources/antismash_db
        antismash --version >> {log}
        antismash --check-prereqs >> {log}
        """

rule antismash:
    input: 
        gbk = "data/processed/genbank/{strains}.gbk",
        resources = "resources/antismash_db/"
    output:
        folder = directory("data/interim/antismash/{version}/{strains}"),
        gbk = "data/interim/antismash/{version}/{strains}/{strains}.gbk",
    conda:
        "../envs/antismash.yaml"
    threads: 8
    log: "workflow/report/{strains}/antismash_{version}.log"
    shell:
        """
        antismash --genefinding-tool prodigal --output-dir {output.folder} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} >> {log}
        """