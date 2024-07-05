rule antismash_db_duckdb_schema:
    input:
        duckdb = "data/processed/{name}/dbt/antiSMASH_{version}/dbt_bgcflow.duckdb"
    output:
        duckdb = temp(directory("data/processed/{name}/antismash_database/antiSMASH_{version}_DB_SCHEMA"))
    conda:
        "../envs/antismash_db-duckdb.yaml"
    log:
        "logs/antismash_db-duckdb_schema/{name}_{version}.log"
    shell:
        """
        duckdb_dir=$PWD/{output.duckdb}
        mkdir -p $duckdb_dir
        duckdb=$PWD/{output.duckdb}/antismash_db.duckdb
        cp {input.duckdb} $duckdb 2>> {log}
        (cd resources/antismash_db-schema_duckdb && python init_duckdb.py --duckdb-database $duckdb --verbose db-schema $duckdb_dir) 2>> {log}
        """

rule antismash_db_duckdb:
    input:
        duckdb = "data/processed/{name}/antismash_database/antiSMASH_{version}_DB_SCHEMA",
        antismash=lambda wildcards: expand("data/processed/{name}/antismash/{version}/{strains}",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
    output:
        database=directory("data/processed/{name}/antismash_database/antiSMASH_database_{version}"),
    conda:
        "../envs/antismash_db-duckdb.yaml"
    threads: 4
    log:
        "logs/antismash_db-duckdb/{name}_{version}.log",
    params:
        prefix=lambda wildcards: find_common_prefix(f"data/processed/{wildcards.name}/antismash/{wildcards.version}/")
    shell:
        """
        set -e
        mkdir -p {output.database}
        bash_script=$PWD/resources/antismash_db-schema_duckdb/full_workflow.sh
        input_dir=$PWD/data/processed/{wildcards.name}/antismash/{wildcards.version}

        # Install the asdb-taxa tool using cargo (Rust's package manager)
        cargo install asdb-taxa &>> {log}

        # Add cargo's bin directory to the PATH to ensure asdb-taxa can be executed
        export PATH="$HOME/.cargo/bin:$PATH"

        # Run the full_workflow.sh script
        duckdb=$PWD/{input.duckdb}/antismash_db.duckdb
        (cd {output.database} && bash $bash_script -p {params.prefix} $input_dir -d $duckdb) &>> {log}
        """
