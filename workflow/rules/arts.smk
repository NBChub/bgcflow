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
        resources="resources/arts",
        khmms="resources/arts/reference/knownresistance.hmm",
    shell:
        """
        mkdir -p data/interim/arts/antismash-{wildcards.version}/{wildcards.strains}
        python {params.resources}/artspipeline1.py {input.antismash} {params.ref} -rd data/interim/arts/antismash-{wildcards.version}/{wildcards.strains} -cpu {threads} -opt kres,phyl,duf -khmms {params.khmms} 2>> {log}
        """

rule arts_extract:
    input:
        arts="data/interim/arts/antismash-{version}/{strains}",
        json="data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        json=temp("data/interim/arts/antismash-{version}/{strains}.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/arts/arts_extract/arts_extract-{version}-{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/arts_extract.py {input.arts}/tables/bgctable.tsv {input.json} {output.json} {wildcards.strains} 2>> {log}
        """

rule arts_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        changelog="data/interim/bgcs/{name}/{version}",
    output:
        table="data/processed/{name}/tables/df_arts_as-{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/arts/arts_combine/arts_combine-{version}-{name}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/arts_gather.py '{input.json}' {input.changelog} {output.table} 2>> {log}
        """
