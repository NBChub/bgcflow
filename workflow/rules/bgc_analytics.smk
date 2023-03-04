rule bgc_count:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.gbk",
    output:
        bgc_count="data/interim/antismash/{version}/{strains}_bgc_counts.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/bgc_counts/as_{version}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_bgc_counts.py {input.antismash} {wildcards.strains} {output.bgc_count} 2>> {log}
        """

rule antismash_overview:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        bgc_table=temp("data/interim/antismash/{version}/{strains}_bgc_overview.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/bgc_overview/as_{version}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_antismash_overview.py {input.antismash} {output.bgc_table} {wildcards.strains} 2>> {log}
        """

rule antismash_overview_gather:
    input:
        bgc_overview=lambda wildcards: expand("data/interim/antismash/{version}/{strains}_bgc_overview.json",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        df_bgc="data/processed/{name}/tables/df_regions_antismash_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/antismash_overview_gather-{version}-{name}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_antismash_overview_gather.py '{input.bgc_overview}' {input.mapping_dir} {output.df_bgc} 2>> {log}
        """

rule antismash_summary:
    input:
        bgc_count=lambda wildcards: expand(
            "data/interim/antismash/{version}/{strains}_bgc_counts.json",
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
        as_dir=lambda wildcards: expand(
            "data/processed/{name}/antismash/{version}/{strains}",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
        bgc_overview="data/processed/{name}/tables/df_regions_antismash_{version}.csv",
    output:
        df_antismash_summary=report(
            "data/processed/{name}/tables/df_antismash_{version}_summary.csv",
            caption="../report/table-antismash.rst",
            category="{name}",
            subcategory="AntiSMASH Summary Table",
        )
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/antismash_summary-{version}-{name}.log",
    params:
        df_samples=lambda wildcards: PEP_PROJECTS[wildcards.name].config["sample_table"],
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_genome_dataset.py '{input.bgc_count}' '{params.df_samples}' {output.df_antismash_summary} 2>> {log}
        """


rule write_dependency_versions:
    output:
        "workflow/report/dependency_versions.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/write_dependency_versions.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_dependencies.py {output} 2> {log}
        """


rule get_project_metadata:
    output:
        "data/processed/{name}/metadata/project_metadata.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/bgc_analytics/get_{name}_metadata.log",
    params:
        bgcflow_version=__version__,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_project_metadata.py {wildcards.name} {output} {params.bgcflow_version} 2> {log}
        """
