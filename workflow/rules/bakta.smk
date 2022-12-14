rule bakta_db_setup:
    output:
        bakta_db = directory("resources/bakta/db"),
    conda:
        "../envs/bakta.yaml"
    log: "workflow/report/logs/bakta/bakta_db_setup/bakta_db_setup.log"
    shell:
        """
        bakta_db download --output /datadrive/bgcflow/resources/bakta
        """

rule bakta:
    input:
        bakta_db = "resources/bakta/db",
        fna = "data/interim/fasta/{strains}.fna",
        org_info = "data/interim/prokka/{strains}/organism_info.txt",
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, "file", DF_SAMPLES, PROKKA_DB_MAP)
    output:
        gff = directory("data/interim/bakta/{strains}"),
        fasta = temp("data/interim/bakta/{strains}.fasta"),
    conda:
        "../envs/bakta.yaml"
    threads: 4
    log: "workflow/report/logs/bakta/bakta/bakta-{strains}.log"
    params:
        gram = "?",
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, "params", DF_SAMPLES, PROKKA_DB_MAP),
    shell:
        """
        cp {input.fna} {output.fasta}
        bakta \
            --output data/interim/bakta/{wildcards.strains} \
            {params.refgbff} \
            --prefix {wildcards.strains} \
            --genus "`cut -d "," -f 1 {input.org_info}`" \
            --species "`cut -d "," -f 2 {input.org_info}`" \
            --strain "`cut -d "," -f 3 {input.org_info}`" \
            --db {input.bakta_db} \
            --gram {params.gram} \
            --threads {threads} \
            --verbose \
            {output.fasta} &> {log}
        """