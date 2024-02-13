rule clone_schema:
    output:
        schema=directory("resources/asdb-schema")
    log:
        "logs/antismash_db/clone_schema.log"
    conda:
        "../envs/antismash.yaml"
    params:
        repository="git@github.com:kblin/asdb-schema.git",
        branch="as7_update"
    shell:
        """
        git clone {params.repository} {output} 2>> {log}
        (cd {output} && git checkout {params.branch}) 2>> {log}
        """

rule download_ncbi_taxdump:
    output:
        ncbi_taxdump="resources/ncbi-taxdump/new_taxdump.tar.gz"
    log:
        "logs/antismash_db/download_ncbi_taxdump.log"
    conda:
        "../envs/antismash.yaml"
    params:
        url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    shell:
        """
        wget -P resources/ncbi-taxdump {params.url} -nc 2>> {log}
        """

rule extract_ncbi_taxdump:
    input:
        rules.download_ncbi_taxdump.output.ncbi_taxdump
    output:
        merged="resources/ncbi-taxdump/merged.dmp",
        rankedlineage="resources/ncbi-taxdump/rankedlineage.dmp"
    log:
        "logs/antismash_db/extract_ncbi_taxdump.log"
    conda:
        "../envs/antismash.yaml"
    shell:
        """
        (cd resources/ncbi-taxdump && tar -xvf new_taxdump.tar.gz) &>> {log}
        """

rule clone_json_importer:
    output:
        import_dir=directory("resources/db-import")
    log:
        "logs/antismash_db/clone_json_importer.log"
    conda:
        "../envs/antismash.yaml"
    params:
        repository="https://github.com/matinnuhamunada/db-import.git",
        branch="atlas-0.1.0-1"
    shell:
        """
        mkdir -p $(dirname {log})
        git clone {params.repository} {output} 2>> {log}
        (cd {output} && git checkout {params.branch}) 2>> {log}
        """

rule asdb_taxa_init:
    input:
        merged = rules.extract_ncbi_taxdump.output.merged,
        rankedlineage = rules.extract_ncbi_taxdump.output.rankedlineage,
        json=lambda wildcards: expand("data/interim/antismash/{version}/{strains}/{strains}.json",
                    name=wildcards.name,
                    version=wildcards.version,
                    strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
                    ),
    output:
        asdb_tax_cache="data/interim/antismash_db/asdb_cache_{name}_{version}.json",
        datadir=temp(directory("data/interim/antismash_db/{name}/{version}"))
    log:
        "logs/antismash_db/asdb_taxa_init_{name}_{version}.log"
    shell:
        """
        mkdir -p {output.datadir}
        for f in {input.json}; do
            cp $f "{output.datadir}/$(basename $f)"
        done
        asdb-taxa init --cache {output.asdb_tax_cache} --datadir {output.datadir} --mergeddump {input.merged} --taxdump {input.rankedlineage} &> {log}
        """

rule import_antismash_db_json:
    input:
        dotenv="config/.env",
        schema=rules.clone_schema.output.schema,
        import_dir=rules.clone_json_importer.output.import_dir,
        taxonomy=rules.asdb_taxa_init.output.asdb_tax_cache,
        json=lambda wildcards: expand("data/interim/antismash/{version}/{strains}/{strains}.json",
                    name=wildcards.name,
                    version=wildcards.version,
                    strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
                    ),
    output:
        db="data/processed/{name}/antismash_db/antismash_db_import_{version}.log"
    log:
        "logs/antismash_db/import_{name}_{version}.log"
    params:
        file_list="data/interim/antismash_db/{name}_{version}_file_list.txt"
    conda:
        "../envs/antismash.yaml"
    shell:
        """
        set -e
        source {input.dotenv}
        (cd {input.schema} && bash init_database.sh) &>> {log}
        echo -e "{input.json}" | tr ' ' '\n' > {params.file_list}
        python {input.import_dir}/import_json.py --taxonomy {input.taxonomy} \
            --from-filelist {params.file_list} \
            --db "host='$PSQL_HOST' port='$PSQL_PORT' user='$PSQL_USER' password='$PGPASSWORD' dbname='$PSQL_DB'" &>> {log}
        echo "Job finished!" > {output.db}
        """
