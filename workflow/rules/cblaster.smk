rule cblaster_genome_db:
    input:
        gbk = lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES, ext="gbk")
    output:
        folder_interim = directory("data/interim/cblaster/{name}/genomes/"),
        folder_processed = directory("data/processed/{name}/cblaster/genomes/"),
        sql = "data/interim/cblaster/{name}/genomes/cblaster_genome_db.sqlite3",
        dmnd = "data/interim/cblaster/{name}/genomes/cblaster_genome_db.dmnd",
        fasta = "data/interim/cblaster/{name}/genomes/cblaster_genome_db.fasta",
    conda:
        "../envs/cblaster.yaml"
    log: "workflow/report/logs/cblaster/cblaster_db_genomes_{name}.log"
    threads: 32
    params:
        db_prefix = "data/interim/cblaster/{name}/genomes/cblaster_genome_db",
        batch_size = 50,
    shell:
        """
        cblaster config --email dummy@cblaster.com
        cblaster makedb --cpus {threads} -b {params.batch_size} -n {params.db_prefix} {input.gbk} 2>> {log}
        cp -r {output.folder_interim} {output.folder_processed} 2>> {log}
        """

rule cblaster_bgc_db:
    input:
        bgc_mapping = "data/interim/bgcs/{name}/{name}_antismash_6.0.1.csv",
    output:
        folder_interim = directory("data/interim/cblaster/{name}/bgcs/"),
        folder_processed = directory("data/processed/{name}/cblaster/bgcs/"),
        sql = "data/interim/cblaster/{name}/bgcs/cblaster_bgc_db.sqlite3",
        dmnd = "data/interim/cblaster/{name}/bgcs/cblaster_bgc_db.dmnd",
        fasta = "data/interim/cblaster/{name}/bgcs/cblaster_bgc_db.fasta", 
    conda:
        "../envs/cblaster.yaml"
    params:
        db_prefix = "data/interim/cblaster/{name}/bgcs/cblaster_bgc_db",
        antismash_dir = "data/interim/bgcs/{name}/6.0.1/*/*region*.gbk",
        batch_size = 50,
    log: "workflow/report/logs/cblaster/cblaster_db_bgc_{name}.log"
    threads: 32
    shell:
        """
        cblaster config --email dummy@cblaster.com
        cblaster makedb --cpus {threads} -b {params.batch_size} -n {params.db_prefix} {params.antismash_dir} 2>> {log}
        cp -r {output.folder_interim} {output.folder_processed} 2>> {log}
        """