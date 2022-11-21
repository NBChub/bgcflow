rule antismash_json_extract:
    input:
        json = "data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        cdss = "data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
        dna_sequences = "data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
        regions = "data/interim/database/as_{version}/{strains}/{strains}_regions.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/scatter/as_{version}_json_extract_{strains}.log"
    params:
        outdir = "data/interim/database/as_{version}/{strains}",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/bgc_meta.py {input.json} {params.outdir} 2>> {log}
        """

rule build_dna_sequences_table:
    input:
        dna_sequences = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_dna_sequences_gather_{name}.log"
    params:
        index_key = "sequence_id"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.dna_sequences}' {params.index_key} {output.dna_sequences} 2>> {log}
        """

rule build_regions_table:
    input:
        regions = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_regions.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_regions_gather_{name}.log"
    params:
        index_key = "region_id"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.regions}' {params.index_key} {output.regions} 2>> {log}
        """

rule build_cdss_table:
    input:
        cdss = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_cdss_gather_{name}.log"
    params:
        index_key = "cds_id"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.cdss}' {params.index_key} {output.cdss} 2>> {log}
        """

### draft rule to build duckdb database
rule build_warehouse:
    input:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        log = "data/processed/{name}/data_warehouse/{version}/database.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/report/database_{version}_{name}.log"
    shell:
        """
        echo {input.cdss} >> {output.log}
        echo {input.regions} >> {output.log}
        echo {input.dna_sequences} >> {output.log}
        """

rule get_dbt_template:
    input:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        #dbt = directory("data/processed/{name}/dbt_as_{version}"),
        profile = "data/processed/{name}/dbt_as_{version}/models/sources.yml"
        # duckdb = "data/processed/{name}/dbt_as_{version}/dbt_bgcflow.duckdb"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "workflow/report/logs/database/report/get_dbt_template_{version}_{name}.log"
    threads: 4
    params:
        dbt = "data/processed/{name}/dbt_as_{version}",
        dbt_repo = "git@github.com:matinnuhamunada/dbt_bgcflow.git",
        cutoff = "0.30",
        as_version = "{version}"
    shell:
        """
        # clone dbt
        if [ -f "{params.dbt}/profiles.yml" ]
        then
            echo "{params.dbt} already exists!" >> {log}
        else
            rm -rf {params.dbt} 2>> {log}
            (cd data/processed/{wildcards.name} \
                && git clone {params.dbt_repo} $(basename {params.dbt})
            ) &>> {log}
        fi

        python {params.dbt}/scripts/source_template.py {params.dbt}/templates/_sources.yml {output.profile} {params.as_version} {params.cutoff} &>> {log}
        """

rule build_database:
    input:
        #dbt = "data/processed/{name}/dbt_as_{version}",
        profile = "data/processed/{name}/dbt_as_{version}/models/sources.yml"
    output:
        duckdb = "data/processed/{name}/dbt_as_{version}/dbt_bgcflow.duckdb"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "workflow/report/logs/database/report/database_{version}_{name}.log"
    threads: 4
    params:
        dbt = "data/processed/{name}/dbt_as_{version}"
    shell:
        """
        (cd {params.dbt} \
            && dbt debug \
            && dbt build \
        ) &>> {log}
        """
