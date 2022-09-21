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