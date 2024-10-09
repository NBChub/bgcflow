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
        if [ {params.release_major} -ge 220 ]; then
            echo "Release major is 220 or above. Using split package option." >> {log}
            # Define the base URL for the split packages
            BASE_URL="https://data.gtdb.ecogenomic.org/releases/release{params.release_major}/{params.release_major}.{params.release_minor}/auxillary_files/gtdbtk_package/split_package"
            # Define the parts
            PARTS=(aa ab ac ad ae af ag ah ai aj ak)
            # Download each part
            for PART in "${{PARTS[@]}}"; do
                wget "$BASE_URL/gtdbtk_r{params.release_major}_data.tar.gz.part_$PART" -O "resources/gtdbtk_r{params.release_major}_data.tar.gz.part_$PART" -nc &>> {log}
            done
            # Concatenate the parts into a single tar.gz file
            cat resources/gtdbtk_r{params.release_major}_data.tar.gz.part_* > resources/gtdbtk_r{params.release_major}_data.tar.gz
            # Extract the concatenated tar.gz file
            mkdir -p resources/gtdbtk
            tar -xvzf resources/gtdbtk_r{params.release_major}_data.tar.gz -C "resources/gtdbtk" --strip 1 &>> {log}
            # Clean up the parts and the concatenated tar.gz file
            rm resources/gtdbtk_r{params.release_major}_data.tar.gz.part_*
            rm resources/gtdbtk_r{params.release_major}_data.tar.gz
        else
            (cd resources && wget https://data.gtdb.ecogenomic.org/releases/release{params.release_major}/{params.release_major}.{params.release_minor}/auxillary_files/gtdbtk_{params.release_version}_data.tar.gz -nc) 2>> {log}
            (cd resources && mkdir -p gtdbtk && tar -xvzf gtdbtk_{params.release_version}_data.tar.gz -C "gtdbtk" --strip 1 && rm gtdbtk_{params.release_version}_data.tar.gz) &>> {log}
        fi
        """

checkpoint prepare_gtdbtk_input:
    input:
        json_list=lambda wildcards: get_json_inputs(wildcards.name, DF_SAMPLES),
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fnadir=directory("data/interim/gtdbtk/{name}/fasta/"),
        fnalist="data/interim/gtdbtk/{name}/fasta_list.txt",
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
        python workflow/bgcflow/bgcflow/data/gtdbtk_prep.py $INPUT_FNA $INPUT_JSON {output.fnadir} {output.fnalist} 2>> {log}
        rm $INPUT_FNA
        rm $INPUT_JSON
        """

rule gtdbtk:
    input:
        gtdbtk="resources/gtdbtk/",
        fnadir="data/interim/gtdbtk/{name}/fasta/",
    output:
        fnalist="data/interim/gtdbtk/{name}/fasta_list_success.txt",
        gtdbtk_dir=directory("data/interim/gtdbtk/{name}/result/"),
        tmpdir=temp(directory("data/interim/gtdbtk/{name}_tmp/")),
        summary_interim="data/interim/gtdbtk/{name}/result/classify/gtdbtk.bac120.summary.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk/gtdbtk/gtdbtk_{name}.log",
    threads: 32
    params:
        ani_screen=ani_screen,
    shell:
        """
        set -e
        mkdir -p {output.tmpdir}
        gtdbtk classify_wf --genome_dir {input.fnadir} --out_dir {output.gtdbtk_dir} --cpus {threads} --pplacer_cpus 1 --tmpdir {output.tmpdir} {params.ani_screen} &>> {log}
        tail -n +2 {output.summary_interim} | cut -f1 > {output.fnalist}
        """

rule gtdbtk_fna_fail:
    output:
        "data/interim/gtdbtk/{name}/fasta_list_fail.txt",
    log:
        "logs/gtdbtk/gtdbtk/gtdbtk_{name}.log",
    shell:
        """
        echo "WARNING: No genomes are eligible for GTDB-Tk classification. Please check if the genome ids already exists in GTDB. Returning empty ouput." > {log}
        echo -e "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings" > {output}
        """

rule evaluate_gtdbtk_input:
    input:
        evaluate_gtdbtk_input
    output:
        summary_processed="data/processed/{name}/tables/gtdbtk.bac120.summary.tsv"
    log:
        "logs/gtdbtk/gtdbtk/evaluate_gtdbtk_{name}.log",
    shell:
        """
        if [ -f "data/interim/gtdbtk/{wildcards.name}/result/classify/gtdbtk.bac120.summary.tsv" ]; then
            cp "data/interim/gtdbtk/{wildcards.name}/result/classify/gtdbtk.bac120.summary.tsv" {output.summary_processed}
        else
            echo "WARNING: No genomes are eligible for GTDB-Tk classification. Please check if the genome ids already exists in GTDB. Returning empty ouput." > {log}
            cp {input} {output.summary_processed}
        fi 2>> {log}
        """
