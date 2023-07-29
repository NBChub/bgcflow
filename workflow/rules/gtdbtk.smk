# Read release version from config
try:
    gtdb_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdb_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdb_release = "207.0"
    gtdb_release_version = "r207_v2"

if "." in str(gtdb_release):
    gtdb_release_major, gtdb_release_minor = str(gtdb_release).split(".")
else:
    gtdb_release_major = gtdb_release
    gtdb_release_minor = "0"

# Decide to use ani screen or not
try:
    if config["rule_parameters"]["gtdbtk"]["ani_screen"]:
        ani_screen = f"--mash_db resources/gtdb-tk-{gtdb_release_version}.msh"
    else:
        ani_screen = "--skip_ani_screen"
except KeyError:
    ani_screen = "--skip_ani_screen"

rule install_gtdbtk:
    output:
        gtdbtk=directory("resources/gtdbtk/"),
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk/gtdbtk-install_gtdbtk.log",
    params:
        release_major=gtdb_release_major,
        release_minor=gtdb_release_minor,
        release_version=gtdb_release_version,
    shell:
        """
        (cd resources && wget https://data.gtdb.ecogenomic.org/releases/release{params.release_major}/{params.release_major}.{params.release_minor}/auxillary_files/gtdbtk_{params.release_version}_data.tar.gz -nc) 2>> {log}
        (cd resources && mkdir -p gtdbtk && tar -xvzf gtdbtk_{params.release_version}_data.tar.gz -C "gtdbtk" --strip 1 && rm gtdbtk_{params.release_version}_data.tar.gz) &>> {log}
        """

rule prepare_gtdbtk_input:
    input:
        json_list=lambda wildcards: get_json_inputs(wildcards.name, DF_SAMPLES),
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fnadir=directory("data/interim/gtdbtk/{name}/fasta/"),
    log:
        "logs/gtdbtk/prepare_gtdbtk_input/{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_FNA="$TMPDIR/df_fna_gtdbtk.txt"
        INPUT_JSON="$TMPDIR/df_json_gtdbtk.txt"
        echo '{input.fna}' > $INPUT_FNA
        echo '{input.json_list}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/gtdbtk_prep.py $INPUT_FNA $INPUT_JSON {output.fnadir} 2>> {log}
        rm $INPUT_FNA
        rm $INPUT_JSON
        """

rule gtdbtk:
    input:
        gtdbtk="resources/gtdbtk/",
        fnadir="data/interim/gtdbtk/{name}/fasta/",
    output:
        #batchfile = "data/interim/gtdbtk/{name}/fasta_batch.tsv",
        gtdbtk_dir=directory("data/interim/gtdbtk/{name}/result/"),
        tmpdir=temp(directory("data/interim/gtdbtk/{name}_tmp/")),
        summary_interim="data/interim/gtdbtk/{name}/result/classify/gtdbtk.bac120.summary.tsv",
        summary_processed="data/processed/{name}/tables/gtdbtk.bac120.summary.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk/gtdbtk/gtdbtk_{name}.log",
    threads: 32
    params:
        ani_screen=ani_screen,
    shell:
        """
        mkdir -p {output.tmpdir}
        gtdbtk classify_wf --genome_dir {input.fnadir} --out_dir {output.gtdbtk_dir} --cpus {threads} --pplacer_cpus 1 --tmpdir {output.tmpdir} {params.ani_screen} &>> {log}
        cp {output.summary_interim} {output.summary_processed}
        """
