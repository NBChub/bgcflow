# Read release version from config
try:
    gtdb_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdb_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdb_release = "207"
    gtdb_release_version = "207_v2"

sys.stderr.write(f"GTDB API | Grabbing metadata using GTDB release version: {gtdb_release_version}\n")

rule gtdb_prep:
    output:
        gtdb_json="data/interim/gtdb/{strains}.json",
    log:
        "logs/gtdb/gtdb_prep/gtdb_prep-{strains}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        samples_path=bgcflow_util_dir / "samples.csv",
        gtdb_paths=GTDB_PATHS,
        version=f"R{str(gtdb_release).split('.')[0]}",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/gtdb_prep.py {wildcards.strains} {output.gtdb_json} '{params.samples_path}' '{params.gtdb_paths}' {params.version} 2> {log}
        """


rule fix_gtdb_taxonomy:
    input:
        json_list=lambda wildcards: [f"data/interim/gtdb/{strains}.json" for strains in PEP_PROJECTS[wildcards.name].sample_tables.genome_id.unique()],
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
        "logs/gtdb/fix_gtdb_taxonomy/fix_gtdb_taxonomy-{name}.log",
    priority: 50
    params:
        samples_path=SAMPLE_PATHS,
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_gtdb.txt"
        echo '{input.json_list}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/fix_gtdb_taxonomy.py $INPUT_JSON {output.meta} 2> {log}
        rm $INPUT_JSON
        """
