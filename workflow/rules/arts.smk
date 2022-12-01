rule arts:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.gbk",
    output:
        folder=directory("data/interim/arts/antismash-{version}/{strains}/"),
    conda:
        "../envs/arts.yaml"
    threads: 4
    log:
        "workflow/report/logs/arts/arts/arts-{version}-{strains}.log",
    params:
        ref="resources/arts/reference/actinobacteria/",
        resources="resources/arts/",
    shell:
        """
        mkdir -p data/interim/arts/antismash-{wildcards.version}/{wildcards.strains}
        python {params.resources}/artspipeline1.py {input.antismash} {params.ref} -rd data/interim/arts/antismash-{wildcards.version}/{wildcards.strains} -cpu {threads} -opt kres,phyl 2>> {log}
        """

rule arts_extract:
    input:
        folder="data/interim/arts/antismash-{version}/{strains}/",
    output:
        json=temp("data/interim/arts/antismash-{version}/{strains}.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/arts/arts/arts_scatter-{version}-{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/arts_extract.py {input.folder}/tables/bgctable.tsv {wildcards.strains} {output.json} 2>> {log}
        """

rule arts_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
    output:
        table="data/processed/{name}/tables/df_arts_as-{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/arts/arts_gather/arts_gather-{version}-{name}.log",
    params:
        index_key = "bgc_id"
    shell:
        """
        python workflow/bgcflow/bgcflow/database/gather_to_csv.py '{input.json}' {params.index_key} {output.table} 2>> {log}
        """
