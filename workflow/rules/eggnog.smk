rule install_eggnog:
    output:
        eggnog_db = directory("resources/eggnog_db"),
        dmnd = "resources/eggnog_db/bacteria.dmnd"
    conda:
        "../envs/eggnog.yaml"
    log: "workflow/report/logs/eggnog-install_eggnog.log"
    shell:
        """
        download_eggnog_data.py --data_dir {output.eggnog_db} -y >> {log}
        create_dbs.py -m diamond --dbname bacteria --taxa Bacteria --data_dir {output.eggnog_db} -y >> {log}
        """

rule eggnog:
    input: 
        faa = "data/interim/prokka/{strains}/{strains}.faa",
        eggnog_db = "resources/eggnog_db",
        dmnd = "resources/eggnog_db/bacteria.dmnd"
    output:
        eggnog_dir = directory("data/interim/eggnog/{strains}/")
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log: "workflow/report/logs/{strains}/eggnog.log"
    shell:
        """
        mkdir -p {output.eggnog_dir}
        emapper.py -i {input.faa} --decorate_gff "yes" --excel --cpu {threads} -o {wildcards.strains} --output_dir {output.eggnog_dir} --data_dir {input.eggnog_db} >> {log}
        """