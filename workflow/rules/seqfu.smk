rule seqfu_stats:
    input: 
        fna = "data/interim/fasta/{strains}.fna"
    output:
        json = "data/interim/seqfu/{strains}.json"
    conda:
        "../envs/seqfu.yaml"
    log: "workflow/report/logs/seqfu/seqfu/seqfu-{strains}.log"
    params:
        precision = 3,
    shell:
        """
        seqfu stats {input.fna} --json -b --gc --precision 3 > {output.json} 2>> {log}
        """

rule seqfu_combine:
    input:
        json = lambda wildcards: expand("data/interim/seqfu/{strains}.json",
                                        strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)])
    output:
        all_csv = report("data/processed/{name}/tables/df_seqfu_stats.csv", caption="../report/table-seqfu.rst", category="{name}", subcategory="Quality Control")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/seqfu/seqfu-{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_seqfu_table.py '{input.json}' {output.all_csv} 2>> {log}
        """

rule seqfu_report:
    input:
        seqfu = "data/processed/{name}/tables/df_seqfu_stats.csv",
        gtdb = "data/processed/{name}/tables/df_gtdb_meta.csv"
    output:
        notebook = "data/processed/{name}/docs/seqfu.ipynb",
        markdown = "data/processed/{name}/docs/seqfu.md",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/seqfu/seqfu-report-{name}.log"
    params:
        notebook = "workflow/notebook/seqfu.py.ipynb"
    shell:
        """
        cp {params.notebook} {output.notebook}
        jupyter nbconvert --to notebook --execute {output.notebook} --output seqfu.ipynb 2>> {log}
        jupyter nbconvert --to markdown {output.notebook} --no-input --output seqfu.md 2>> {log}
        """
