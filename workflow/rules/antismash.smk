rule antismash_db_setup:
    output:
        directory("resources/antismash_db"),
    conda:
        "../envs/antismash.yaml"
    log:
        "workflow/report/logs/antismash/antismash_db_setup.log",
    shell:
        """
        download-antismash-databases --database-dir resources/antismash_db
        antismash --version >> {log}
        antismash --check-prereqs >> {log}
        """


rule antismash:
    input:
        gbk="data/interim/processed-genbank/{strains}.gbk",
        resources="resources/antismash_db/",
    output:
        folder=directory("data/interim/antismash/{version}/{strains}/"),
        gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
        json="data/interim/antismash/{version}/{strains}/{strains}.json",
        zip="data/interim/antismash/{version}/{strains}/{strains}.zip",
    conda:
        "../envs/antismash.yaml"
    threads: 4
    log:
        "workflow/report/logs/antismash/antismash/antismash_{version}-{strains}.log",
    params:
        folder=directory("data/interim/antismash/{version}/{strains}/"),
        genefinding="none",
    shell:
        """
        antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} --logfile {log} 2>> {log}
        """


rule copy_antismash:
    input:
        dir="data/interim/antismash/{version}/{strains}",
    output:
        dir=directory("data/processed/{name}/antismash/{version}/{strains}"),
    conda:
        "../envs/antismash.yaml"
    log:
        "workflow/report/logs/antismash/copy_antismash/copy_antismash_{version}-{strains}-{name}.log",
    shell:
        """
        base_dir=$PWD
        mkdir {output.dir}
        (cd {output.dir} && for item in $(ls $base_dir/{input.dir}); do ln -s $base_dir/{input.dir}/$item $(basename $item); done) 2>> {log}
        """
