rule arts:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.gbk",
    output:
        folder=directory("data/interim/arts/antismash-{version}/{strains}/"),
    conda:
        "../envs/arts.yaml"
    threads: 4
    log:
        "logs/arts/arts/arts-{version}-{strains}.log",
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
        all_hits=temp("data/interim/arts/antismash-{version}/{strains}_arts_all_hits.json"),
        bgctable=temp("data/interim/arts/antismash-{version}/{strains}_arts_bgctable_summary.json"),
        coretable=temp("data/interim/arts/antismash-{version}/{strains}_arts_coretable_summary.json"),
        knownhits=temp("data/interim/arts/antismash-{version}/{strains}_arts_knownhits.json")
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_extract/arts_extract-{version}-{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/arts_extract_all.py {input.arts} {wildcards.strains} {input.json} data/interim/arts/antismash-{wildcards.version}  2>> {log}
        """

rule arts_allhits_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}_arts_all_hits.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        changelog="data/interim/bgcs/{name}/{version}",
    output:
        table=temp("data/processed/{name}/tables/df_arts_all_hits_as-{version}.csv")
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_combine/arts_combine_all_hits-{version}-{name}.log",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_arts_allhits.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/arts_gather.py $INPUT_JSON {input.changelog} {output.table} 2>> {log}
        rm $INPUT_JSON
        """

rule arts_bgctable_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}_arts_bgctable_summary.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        changelog="data/interim/bgcs/{name}/{version}",
    output:
        table="data/processed/{name}/tables/df_arts_bgctable_as-{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_combine/arts_combine_bgctable-{version}-{name}.log",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_arts_bgctable.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/arts_gather.py $INPUT_JSON {input.changelog} {output.table} 2>> {log}
        rm $INPUT_JSON
        """

rule arts_coretable_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}_arts_coretable_summary.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        changelog="data/interim/bgcs/{name}/{version}",
    output:
        table="data/processed/{name}/tables/df_arts_coretable_as-{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_combine/arts_combine_coretable-{version}-{name}.log",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_arts_coretable.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/arts_gather.py $INPUT_JSON {input.changelog} {output.table} 2>> {log}
        rm $INPUT_JSON
        """

rule arts_knownhits_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/arts/antismash-{version}/{strains}_arts_knownhits.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
            version=wildcards.version
        ),
        changelog="data/interim/bgcs/{name}/{version}",
    output:
        table="data/processed/{name}/tables/df_arts_knownhits_as-{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_combine/arts_combine_knownhits-{version}-{name}.log",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_arts_knownhits.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/arts_gather.py $INPUT_JSON {input.changelog} {output.table} 2>> {log}
        rm $INPUT_JSON
        """

rule arts_final:
    input:
        all_hits="data/processed/{name}/tables/df_arts_all_hits_as-{version}.csv",
        bgctable="data/processed/{name}/tables/df_arts_bgctable_as-{version}.csv",
        coretable="data/processed/{name}/tables/df_arts_coretable_as-{version}.csv",
        knownhits="data/processed/{name}/tables/df_arts_knownhits_as-{version}.csv"
    output:
        table="data/processed/{name}/tables/df_arts_as-{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/arts/arts_combine/arts_combine_final-{version}-{name}.log",
    shell:
        """
        cp {input.all_hits} {output.table} 2>> {log}
        """