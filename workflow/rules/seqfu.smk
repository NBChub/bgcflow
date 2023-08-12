rule seqfu_stats:
    input:
        fna="data/interim/fasta/{strains}.fna",
    output:
        json="data/interim/seqfu/{strains}.json",
    conda:
        "../envs/seqfu.yaml"
    log:
        "logs/seqfu/seqfu/seqfu-{strains}.log",
    params:
        precision=3,
    shell:
        """
        seqfu stats {input.fna} --json -b --gc --precision 3 > {output.json} 2>> {log}
        """


rule seqfu_combine:
    input:
        json=lambda wildcards: expand(
            "data/interim/seqfu/{strains}.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        all_csv=report(
            "data/processed/{name}/tables/df_seqfu_stats.csv",
            caption="../report/table-seqfu.rst",
            category="{name}",
            subcategory="Quality Control",
        ),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/seqfu/seqfu-{name}.log",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_seqfu.txt"
        echo '{input.json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/make_seqfu_table.py $INPUT_JSON {output.all_csv} &>> {log}
        rm $INPUT_JSON
        """
