# Read release version from config
try:
    gtdbtk_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdbtk_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdbtk_release = "207"
    gtdbtk_release_version = "207_v2"

sys.stderr.write(f"Using GTDB-tk release version: {gtdbtk_release_version}\n")


rule install_gtdbtk:
    output:
        gtdbtk=directory("resources/gtdbtk/"),
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "workflow/report/logs/gtdbtk/gtdbtk-install_gtdbtk.log",
    params:
        release=gtdbtk_release,
        release_version=gtdbtk_release_version,
    shell:
        """
        (cd resources && wget https://data.gtdb.ecogenomic.org/releases/release{params.release}/{params.release}.0/auxillary_files/gtdbtk_r{params.release_version}_data.tar.gz -nc) 2>> {log}
        (cd resources && mkdir -p gtdbtk && tar -xvzf gtdbtk_r{params.release_version}_data.tar.gz -C "gtdbtk" --strip 1 && rm gtdbtk_r{params.release_version}_data.tar.gz) 2>> {log}
        """


rule prepare_gtdbtk_input:
    input:
        json_list=lambda wildcards: get_json_inputs(wildcards.name, DF_SAMPLES),
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fnadir=directory("data/interim/gtdbtk/{name}/fasta/"),
    log:
        "workflow/report/logs/gtdbtk/prepare_gtdbtk_input/{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/gtdbtk_prep.py '{input.fna}' '{input.json_list}' {output.fnadir} 2>> {log}
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
        "workflow/report/logs/gtdbtk/gtdbtk/gtdbtk_{name}.log",
    threads: 32
    shell:
        """
        mkdir -p {output.tmpdir}
        gtdbtk classify_wf --genome_dir {input.fnadir} --out_dir {output.gtdbtk_dir} --cpus {threads} --pplacer_cpus 1 --tmpdir {output.tmpdir} &>> {log}
        cp {output.summary_interim} {output.summary_processed}
        """
