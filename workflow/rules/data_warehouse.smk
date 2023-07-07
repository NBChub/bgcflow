rule csv_to_parquet:
    input:
        csv=final_outputs,
    output:
        parquet="data/processed/{name}/data_warehouse/tables.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/data_warehouse/{name}/convert_to_parquet.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/database/csv_to_parquet.py data/processed/{wildcards.name} 2>> {log}
        touch {output.parquet}
        """
