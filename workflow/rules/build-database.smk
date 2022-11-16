rule antismash_json_extract:
    input:
        json="data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        cdss="data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
        dna_sequences="data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
        regions="data/interim/database/as_{version}/{strains}/{strains}_regions.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/database/scatter/as_{version}_json_extract_{strains}.log",
    params:
        outdir="data/interim/database/as_{version}/{strains}",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/bgc_meta.py {input.json} {params.outdir} 2>> {log}
        """


rule build_dna_sequences_table:
    input:
        dna_sequences=lambda wildcards: expand(
            "data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        dna_sequences="data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/database/gather/as_{version}_dna_sequences_gather_{name}.log",
    params:
        index_key="sequence_id",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.dna_sequences}' {params.index_key} {output.dna_sequences} 2>> {log}
        """


rule build_regions_table:
    input:
        regions=lambda wildcards: expand(
            "data/interim/database/as_{version}/{strains}/{strains}_regions.json",
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        regions="data/processed/{name}/data_warehouse/{version}/regions.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/database/gather/as_{version}_regions_gather_{name}.log",
    params:
        index_key="region_id",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.regions}' {params.index_key} {output.regions} 2>> {log}
        """


rule build_cdss_table:
    input:
        cdss=lambda wildcards: expand(
            "data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        cdss="data/processed/{name}/data_warehouse/{version}/cdss.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/database/gather/as_{version}_cdss_gather_{name}.log",
    params:
        index_key="cds_id",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.cdss}' {params.index_key} {output.cdss} 2>> {log}
        """


### draft rule to build duckdb database
rule build_warehouse:
    input:
        cdss="data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions="data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences="data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        log="data/processed/{name}/data_warehouse/{version}/database.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/database/report/database_{version}_{name}.log",
    shell:
        """
        echo {input.cdss} >> {output.log}
        echo {input.regions} >> {output.log}
        echo {input.dna_sequences} >> {output.log}
        """


rule clone_dbt:
    output:
        dbt=directory("resources/bgcflow_dbt-duckdb"),
    conda:
        "../envs/dbt-duckdb.yaml"
    log:
        "workflow/report/logs/database/clone_dbt-duckdb.log",
    params:
        dbt_repo="git@github.com:NBChub/bgcflow_dbt-duckdb.git",
    shell:
        """
        # clone dbt
        if [ -f "{output.dbt}/profiles.yml" ]
        then
            echo "{output.dbt} already exists!" >> {log}
        else
            rm -rf {output.dbt} 2>> {log}
            (cd resources \
                && git clone {params.dbt_repo} $(basename {output.dbt})
            ) &>> {log}
        fi
        """


rule build_database:
    input:
        repo="resources/bgcflow_dbt-duckdb",
        cdss="data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions="data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences="data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        dbt=directory("data/processed/{name}/dbt_as_{version}"),
        duckdb="data/processed/{name}/dbt_as_{version}/dbt_bgcflow.duckdb",
    conda:
        "../envs/dbt-duckdb.yaml"
    log:
        "workflow/report/logs/database/report/database_{version}_{name}.log",
    threads: 4
    shell:
        """
        cp -r {input.repo}/* {output.dbt}/. 2>> {log}

        # run dbt
        (cd {output.dbt} \
            && dbt debug --profiles-dir . \
            && dbt build --profiles-dir .
        ) &>> {log}
        """
