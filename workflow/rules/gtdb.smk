import requests
import multiprocessing
import time

# Read release version from config
try:
    gtdb_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdb_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdb_release = "207"
    gtdb_release_version = "207_v2"
sys.stderr.write(f"Checking GTDB API...\n")
sys.stderr.write(f" - GTDB API | Grabbing metadata using GTDB release version: {gtdb_release_version}\n")

# Testing connection to GTDB API
gtdb_api_url = "https://gtdb-api.ecogenomic.org/status/db"

gtdb_offline_mode = False  # Flag to indicate if the GTDB API is offline

if "rule_parameters" in config.keys():
    if "use_gtdb_api" in config["rule_parameters"]:
        if config["rule_parameters"]["use_gtdb_api"] == False:
            gtdb_offline_mode = True

if not gtdb_offline_mode:
    try:
        sys.stderr.write(f" - GTDB API | Testing connection to: {gtdb_api_url}\n")
        response = requests.get(gtdb_api_url)
        response.raise_for_status()  # Raise an exception for HTTP errors (4xx, 5xx)
        data = response.json()  # Assuming the API returns JSON data
        sys.stderr.write(f' - GTDB API | Database is online: {data["online"]}\n')

    except requests.exceptions.RequestException as e:
        sys.stderr.write(f" - GTDB API | Error: {e}\n")
        sys.stderr.write(" - GTDB API | Failed to connect to the GTDB API.\n")
        sys.stderr.write(" - GTDB API | It is possible to continue in offline mode. This will return empty taxonomic information for all NCBI genomes!\n")

        # Check if the error is due to a 504 or 404 status code
        if response.status_code == 504 or response.status_code == 404:
            countdown_seconds = 30
            while countdown_seconds > 0:
                sys.stderr.write(f"\rDo you want to continue in offline mode? (yes/no/stop) (Time left: {countdown_seconds} seconds): ")
                sys.stderr.flush()
                time.sleep(1)  # Wait for 1 second
                countdown_seconds -= 1

            sys.stderr.write("\rDo you want to continue in offline mode? (yes/no/stop) (Time left: 0 seconds): \n")
            sys.stderr.flush()

            user_input = get_user_input_with_timeout("", 0)  # Get user input with no timeout
            if user_input is not None and user_input.strip().lower() == 'yes':
                gtdb_offline_mode = True  # Continue in offline mode
            elif user_input is not None and user_input.strip().lower() == 'stop':
                raise
            else:
                # timed out, continue anyway
                #sys.stderr.write("WARNING: No response, continuing BGCFlow in online mode anyway...")
                pass

else:
    sys.stderr.write("WARNING! GTDB API offline mode enabled. This will return empty taxonomic information for all NCBI genomes!\n")

sys.stderr.write(f" - GTDB API | Searching in offline mode: {gtdb_offline_mode}\n")

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
        gtdb_link = f"https://data.gtdb.ecogenomic.org/releases/release{str(gtdb_release).split('.')[0]}/{gtdb_release}/bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv.gz",
        gtdb_table = f"bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv",
        offline=gtdb_offline_mode,
        api_base="https://gtdb-api.ecogenomic.org",
    shell:
        """
        if python workflow/bgcflow/bgcflow/data/gtdb_prep.py {wildcards.strains} {output.gtdb_json} '{params.samples_path}' '{params.gtdb_paths}' {params.version} {params.offline} {params.api_base} 2> {log}; then
            echo "gtdb_prep.py executed successfully" >> {log}
        else
            echo "gtdb_prep.py failed, getting dataset from table instead..." >> {log}
            if [ ! -f resources/gtdb_download/{params.gtdb_table} ]; then
                echo "Downloading GTDB table from {params.gtdb_link}" >> {log}
                mkdir -p resources/gtdb_download/
                wget -P resources/gtdb_download/ {params.gtdb_link} -nc 2>> {log}
                gunzip -c resources/gtdb_download/{params.gtdb_table}.gz > resources/gtdb_download/{params.gtdb_table} 2>> {log}
            fi
            echo "Running gtdb_prep_from_table.py" >> {log}
            python workflow/bgcflow/bgcflow/data/gtdb_prep_from_table.py {wildcards.strains} resources/gtdb_download/{params.gtdb_table} {params.version} {output.gtdb_json} 2>> {log}
        fi
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
