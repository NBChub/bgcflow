# Read release version from config
try:
    gtdbtk_release = config['rule_parameters']['install_gtdbtk']['release']
    gtdbtk_release_version = config['rule_parameters']['install_gtdbtk']['release_version']
except KeyError:
    gtdbtk_release = '220'
    gtdbtk_release_version = '220'

sys.stderr.write(f'Using GTDB-tk release version: {gtdbtk_release_version}\n')

rule install_gtdbtk:
    output:
        checkm_db = directory("resources/gtdbtk/")
    conda:
        "../envs/gtdbtk.yaml"
    log: "workflow/report/logs/gtdbtk/gtdbtk-install_gtdbtk.log"
    params:
        release = gtdbtk_release,
        release_version = gtdbtk_release_version
    shell:
        """
        (cd resources && wget https://data.gtdb.ecogenomic.org/releases/release{params.release}/{params.release}.0/auxillary_files/gtdbtk_r{params.release_version}_data.tar.gz) 2>> {log}
        (cd resources && mkdir -p gtdbtk && tar -xvzf gtdbtk_r{params.release_version}_data.tar.gz -C "gtdbtk" --strip 1 && rm gtdbtk_{params.release_version}_data.tar.gz) 2>> {log}
        conda env config vars set GTDBTK_DATA_PATH="resources/gtdbtk" 2>> {log}
        """

rule gtdbtk:
    input:
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
        checkm_db = "resources/gtdbtk/",
    output:
        batchfile = "data/interim/gtdbtk/{name}/fasta_batch.tsv",
        gtdbtk_dir = directory("data/interim/gtdbtk/{name}/"),
        tmpdir = temp(directory("data/interim/gtdbtk/{name}_tmp/")),
        fnadir = temp(directory("data/interim/gtdbtk/{name}/fasta/")),
        summary_interim = "data/interim/gtdbtk/{name}/classify/gtdbtk.bac120.summary.tsv",
        summary_processed = "data/processed/{name}/tables/gtdbtk.bac120.summary.tsv",
    conda:
        "../envs/gtdbtk.yaml"
    log: "workflow/report/logs/gtdbtk/gtdbtk_{name}.log"
    threads: 32
    shell:
        """
        mkdir -p {output.fnadir}
        for fna in {input.fna} 
        do
            cp {input.fna} {output.fnadir} 
        done
        mkdir -p {output.tmpdir}
        gtdbtk classify_wf --genome_dir {output.fnadir} --out_dir {output.gtdbtk_dir} --cpus {threads} --pplacer_cpus 1 --tmpdir {output.tmpdir}
        cp {output.summary_interim} {output.summary_processed}
        """
