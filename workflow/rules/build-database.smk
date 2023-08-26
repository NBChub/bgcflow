rule antismash_json_extract:
    input:
        json = "data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        cdss = temp("data/interim/database/as_{version}/{strains}/{strains}_cdss.json"),
        dna_sequences = temp("data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json"),
        regions = temp("data/interim/database/as_{version}/{strains}/{strains}_regions.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/scatter/as_{version}_json_extract_{strains}.log"
    params:
        outdir = "data/interim/database/as_{version}/{strains}",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/bgc_meta.py {input.json} {params.outdir} {wildcards.strains} 2>> {log}
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
    log: "logs/database/gather/as_{version}_dna_sequences_gather_{name}.log"
    params:
        index_key = "sequence_id",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet.py '{input.dna_sequences}' {params.index_key} {output.dna_sequences} 2>> {log}
        """

rule build_regions_table:
    input:
        regions = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_regions.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/gather/as_{version}_regions_gather_{name}.log"
    params:
        index_key = "region_id",
        exclude = "_regions.json"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet_with_correction.py '{input.regions}' {input.mapping_dir} {output.regions} {params.exclude} {params.index_key} 2>> {log}
 2>> {log}
        """

rule build_cdss_table:
    input:
        cdss = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/gather/as_{version}_cdss_gather_{name}.log"
    params:
        index_key = "cds_id",
        exclude = "_cdss.json"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_parquet_with_correction.py '{input.cdss}' {input.mapping_dir} {output.cdss} {params.exclude} {params.index_key} 2>> {log}
 2>> {log}
        """

rule build_warehouse:
    input:
        cdss = "data/processed/{name}/data_warehouse/{version}/cdss.parquet",
        regions = "data/processed/{name}/data_warehouse/{version}/regions.parquet",
        dna_sequences = "data/processed/{name}/data_warehouse/{version}/dna_sequences.parquet",
    output:
        log = "data/processed/{name}/data_warehouse/{version}/database.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/database/report/database_{version}_{name}.log"
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
        profile = "data/processed/{name}/dbt/antiSMASH_{version}/models/sources.yml"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "logs/database/report/get_dbt_template_{version}_{name}.log"
    threads: 4
    params:
        dbt = "data/processed/{name}/dbt/antiSMASH_{version}",
        dbt_repo = "https://github.com/NBChub/bgcflow_dbt-duckdb",
        release = "0.2.1",
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
            mkdir -p data/processed/{wildcards.name}/dbt
            (cd data/processed/{wildcards.name}/dbt \
                && wget {params.dbt_repo}/archive/refs/tags/v{params.release}.zip && unzip v{params.release}.zip
            ) &>> {log}
            mv data/processed/{wildcards.name}/dbt/bgcflow_dbt-duckdb-{params.release} {params.dbt}
            rm data/processed/{wildcards.name}/dbt/v{params.release}.zip
        fi

        python {params.dbt}/scripts/source_template.py {params.dbt}/templates/_sources.yml {output.profile} {params.as_version} {params.cutoff} &>> {log}
        """

def exclude_model_dbt(model_to_ignore):
    """
    Returns a string containing the `--exclude` option followed by the models to ignore in a dbt project.

    Args:
        model_to_ignore (list of str): A list of model names to ignore.

    Returns:
        str: A string containing the `--exclude` option followed by the model names, or an empty string if the list is empty.

    Example:
        >>> exclude_model_dbt(['model1', 'model2'])
        '--exclude model1 model2'

    Note:
        The returned string can be used as an argument to the `dbt build` command to exclude the specified models from the build process.
    """
    if len(model_to_ignore) == 0:
        return ""
    else:
        return " ".join(["--exclude"] + model_to_ignore)

rule build_database:
    input:
        profile = "data/processed/{name}/dbt/antiSMASH_{version}/models/sources.yml"
    output:
        duckdb = "data/processed/{name}/dbt/antiSMASH_{version}/dbt_bgcflow.duckdb"
    conda:
        "../envs/dbt-duckdb.yaml"
    log: "logs/database/report/database_{version}_{name}.log"
    threads: 16
    params:
        dbt = "data/processed/{name}/dbt/antiSMASH_{version}",
        exclude = lambda wildcards: exclude_model_dbt(models_to_ignore[wildcards.name])
    shell:
        """
        command="dbt build --threads {threads} {params.exclude} -x"
        echo $command >> {log}
        (cd {params.dbt} \
            && dbt debug \
            && $command \
        ) &>> {log}
        """
