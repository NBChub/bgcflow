rule metabase_install:
    output:
        jar="resources/metabase/metabase.jar",
    log:
        "logs/metabase_install.log",
    conda:
        "../envs/utilities.yaml"
    params:
        version="v0.45.3"
    shell:
        """
        wget -O {output.jar} https://downloads.metabase.com/{params.version}/metabase.jar 2>> {log}
        """

rule metabase_duckdb_plugin:
    output:
        plugin="resources/metabase/plugins/duckdb.metabase-driver.jar",
    log:
        "logs/metabase_duckdb_install.log",
    params:
        release="0.1.6"
    shell:
        """
        wget -O {output.plugin} https://github.com/AlexR2D2/metabase_duckdb_driver/releases/download/{params.release}/duckdb.metabase-driver.jar 2>> {log}
        """
