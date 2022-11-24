# Read release version from config
try:
    gtdbtk_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdbtk_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdbtk_release = "207"
    gtdbtk_release_version = "207_v2"


rule gtdb_prep:
    output:
        gtdb_json="data/interim/gtdb/{strains}.json",
    log:
        "workflow/report/logs/gtdb/gtdb_prep/gtdb_prep-{strains}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        samples_path=bgcflow_util_dir / "samples.csv",
        gtdb_paths=GTDB_PATHS,
        version=f"R{gtdbtk_release}",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/gtdb_prep.py {wildcards.strains} {output.gtdb_json} '{params.samples_path}' '{params.gtdb_paths}' {params.version} 2> {log}
        """


rule fix_gtdb_taxonomy:
    input:
        json_list=lambda wildcards: get_json_inputs(wildcards.name, DF_SAMPLES),
    output:
        meta=report(
            "data/processed/{name}/tables/df_gtdb_meta.csv",
            caption="../report/table-gtdb.rst",
            category="{name}",
            subcategory="Taxonomic Placement",
        ),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "workflow/report/logs/gtdb/fix_gtdb_taxonomy/fix_gtdb_taxonomy-{name}.log",
    priority: 50
    params:
        samples_path=SAMPLE_PATHS,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/fix_gtdb_taxonomy.py '{input.json_list}' {output.meta} 2> {log}
        """
