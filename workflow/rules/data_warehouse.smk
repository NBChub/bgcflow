rule copy_mibig_table:
    input:
        "resources/mibig/df_mibig_bgcs.csv"
    output:
        temp("data/processed/{name}/tables/df_mibig_bgcs.csv")
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/data_warehouse/{name}/copy_mibig.log",
    shell:
        """
        cp {input} {output} 2>> {log}
        """

rule csv_to_parquet:
    input:
        csv=final_outputs,
        mibig="data/processed/{name}/tables/df_mibig_bgcs.csv"
    output:
        parquet="data/processed/{name}/data_warehouse/tables/df_mibig_bgcs.parquet",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/data_warehouse/{name}/convert_to_parquet.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/csv_to_parquet.py data/processed/{wildcards.name} 2>> {log}
        touch {output.parquet}
        """
