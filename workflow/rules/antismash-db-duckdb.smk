rule antismash_db_duckdb:
    input:
        antismash=lambda wildcards: expand("data/processed/{name}/antismash/{version}/{strains}",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
    output:
        database=directory("data/processed/{name}/dbt/antiSMASH_database_{version}"),
    conda:
        "../envs/antismash_db-duckdb.yaml"
    threads: 4
    log:
        "logs/antismash_db-duckdb/{name}_{version}.log",
    shell:
        """
        PWD=$(pwd)
        mkdir -p {output.database}
        bash_script=$PWD/resources/antismash_db-schema_duckdb/full_workflow.sh
        input_dir=$PWD/data/processed/{wildcards.name}/antismash/{wildcards.version}

        # Install the asdb-taxa tool using cargo (Rust's package manager)
        cargo install asdb-taxa 2>> {log}

        # Add cargo's bin directory to the PATH to ensure asdb-taxa can be executed
        export PATH="$HOME/.cargo/bin:$PATH"

        # Run the full_workflow.sh script
        (cd {output.database} && bash $bash_script $input_dir) &>> {log}
        """
