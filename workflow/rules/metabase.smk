rule metabase_install:
    output:
        jar="resources/metabase/metabase_{METABASE_VERSION}.jar",
    log:
        "logs/metabase/metabase_install_{METABASE_VERSION}.log",
    conda:
        "../envs/utilities.yaml"
    params:
        version=metabase_config["METABASE_VERSION"]
    shell:
        """
        wget -O {output.jar} https://downloads.metabase.com/{params.version}/metabase.jar 2>> {log}
        """

rule metabase_duckdb_plugin:
    output:
        plugin="resources/metabase/plugins/duckdb.metabase-driver_{METABASE_DUCKDB_PLUGIN_VERSION}.jar",
    log:
        "logs/metabase/metabase_duckdb_install_{METABASE_DUCKDB_PLUGIN_VERSION}.log",
    params:
        release=metabase_config["METABASE_DUCKDB_PLUGIN_VERSION"]
    shell:
        """
        wget -O {output.plugin} https://github.com/AlexR2D2/metabase_duckdb_driver/releases/download/{params.release}/duckdb.metabase-driver.jar 2>> {log}
        """
