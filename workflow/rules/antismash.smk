rule antismash_db_setup:
    output:
        directory("resources/antismash_db"),
    conda:
        "../envs/antismash.yaml"
    log: "workflow/report/logs/antismash_db_setup.log"
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
        gbk = "data/interim/antismash/{version}/{strains}/{strains}.gbk",
    conda:
        "../envs/antismash.yaml"
    threads: 8
    log: "workflow/report/logs/{strains}/antismash_{version}.log"
    params:
        folder = directory("data/interim/antismash/{version}/{strains}"),
    shell:
        """
        antismash --genefinding-tool prodigal --output-dir {params.folder} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} -v
        ls {params.folder} > {log}
        """