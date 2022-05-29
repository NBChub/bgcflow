rule install_checkm:
    output:
        checkm_db = directory("resources/checkm/")
    conda:
        "../envs/checkm.yaml"
    log: "workflow/report/logs/checkm/checkm-install_checkm.log"
    shell:
        """
        (cd resources && wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz) 2>> {log}
        (cd resources && mkdir -p checkm && tar -xvf checkm_data_2015_01_16.tar.gz -C checkm && rm checkm_data_2015_01_16.tar.gz) 2>> {log}
        checkm data setRoot resources/checkm 2>> {log}
        """

rule checkm:
    input:
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
        checkm_db = "resources/checkm/",
    output:
        fna = temp(directory('data/interim/checkm/{name}_fna')),
        stat = "data/interim/checkm/{name}/storage/bin_stats_ext.tsv",
        checkm_dir = directory("data/interim/checkm/{name}"),
    conda:
        "../envs/checkm.yaml"
    log: "workflow/report/logs/checkm/checkm_{name}.log"
    params: 
        checkm_log = "data/interim/checkm/{name}/checkm_{name}.log",
    threads: 16
    shell:
        """
        mkdir -p {output.fna}
        for f in {input.fna}; do cp $f {output.fna}/.; done
        checkm lineage_wf -t {threads} --reduced_tree -x fna {output.fna} {output.checkm_dir} &>> {log}
        """

rule checkm_out:
    input:
        stat = "data/interim/checkm/{name}/storage/bin_stats_ext.tsv",
    output:
        stat_processed = report("data/processed/{name}/tables/df_checkm_stats.csv",  caption="../report/table-checkm.rst", category="Quality Control"),
    log: "workflow/report/logs/checkm/checkm_out_{name}.log"
    params:
        checkm_json = directory("data/interim/checkm/{name}/json/"),
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        mkdir -p {params.checkm_json}
        python workflow/bgcflow/bgcflow/data/get_checkm_data.py {input.stat} {params.checkm_json} {output.stat_processed} 2>> {log}
        """