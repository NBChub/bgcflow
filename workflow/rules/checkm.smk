rule install_checkm:
    output:
        checkm_db = directory("resources/checkm/")
    conda:
        "../envs/checkm.yaml"
    log: "workflow/report/logs/checkm/checkm-install_checkm.log"
    shell:
        """
        (cd resources && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz) 2> {log}
        (cd resources && mkdir -p checkm && tar -xvf checkm_data_2015_01_16.tar.gz -C checkm && rm checkm_data_2015_01_16.tar.gz) 2>> {log}
        checkm data setRoot resources/checkm >> {log}
        """

rule checkm:
    input:
        _all_ = expand("data/interim/fasta/{strains}.fna", strains=STRAINS),
        fna = "data/interim/fasta",
        checkm_db = "resources/checkm/",
    output:
        stat = "data/interim/checkm/storage/bin_stats_ext.tsv",
        checkm_dir = directory("data/interim/checkm/"),
        bins = expand("data/interim/checkm/bins/{strains}/genes.faa", strains=STRAINS),
    conda:
        "../envs/checkm.yaml"
    log: "workflow/report/logs/checkm/checkm.log"
    params: 
        checkm_log = "data/interim/checkm/checkm.log",
    threads: 16
    shell:
        """
        checkm lineage_wf -t {threads} --reduced_tree -x fna {input.fna} {output.checkm_dir} 2>> {log}
        cat {params.checkm_log} >> {log}
        """

rule checkm_out:
    input:
        stat = "data/interim/checkm/storage/bin_stats_ext.tsv",
        _all_ = expand("data/interim/checkm/bins/{strains}/genes.faa", strains=STRAINS),
    output:
        stat_processed = report("data/processed/tables/df_checkm_stats.csv",  caption="../report/table-checkm.rst", category="Quality Control"),
        checkm_json = directory("data/interim/checkm/json/"),
        json_all = expand("data/interim/checkm/json/{strains}.json", strains=STRAINS),
    log: "workflow/report/logs/checkm/checkm_out.log"
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_checkm_data.py {input.stat} {output.checkm_json} {output.stat_processed} 2>> {log}
        """