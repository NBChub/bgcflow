def find_common_prefix(input_dir):
    """
    Find common prefixes of antiSMASH JSON files in the input directory.
    """
    input_path = Path(input_dir)

    # Find all .json files in the input directory recursively
    json_files = list(input_path.glob('**/*.json'))

    # Group files by the first letter of each file
    files_by_first_letter = defaultdict(list)
    for file in json_files:
        first_letter = file.stem[0] if file.stem else ''
        files_by_first_letter[first_letter].append(file)

    # Find the common prefix within each group
    common_prefixes = []
    for files in files_by_first_letter.values():
        filenames = [file.stem for file in files]
        common_prefix = os.path.commonprefix(filenames)
        if common_prefix:  # Only add non-empty prefixes
            common_prefixes.append(common_prefix)

    # Join the common prefixes with ","
    return ",".join(common_prefixes)

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
    params:
        prefix=lambda wildcards: find_common_prefix(f"data/processed/{wildcards.name}/antismash/{wildcards.version}/")
    shell:
        """
        set -e
        PWD=$(pwd)
        mkdir -p {output.database}
        bash_script=$PWD/resources/antismash_db-schema_duckdb/full_workflow.sh
        input_dir=$PWD/data/processed/{wildcards.name}/antismash/{wildcards.version}

        # Install the asdb-taxa tool using cargo (Rust's package manager)
        cargo install asdb-taxa &>> {log}

        # Add cargo's bin directory to the PATH to ensure asdb-taxa can be executed
        export PATH="$HOME/.cargo/bin:$PATH"

        # Run the full_workflow.sh script
        (cd {output.database} && bash $bash_script -p {params.prefix} $input_dir) &>> {log}
        """
