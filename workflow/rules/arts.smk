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

rule arts_summarize:
    input:
        folder=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
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
        python workflow/bgcflow/bgcflow/data/arts_extract.py '{input.folder}' {input.bgc_mapping} {output.table} 2>> {log}
        """
