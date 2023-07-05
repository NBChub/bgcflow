antismash_major_version = int(dependency_version["antismash"].split(".")[0])

if antismash_major_version <= 6:
    rule antismash_db_setup:
        output:
            directory("resources/antismash_db"),
        conda:
            "../envs/antismash.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        shell:
            """
            download-antismash-databases --database-dir resources/antismash_db 2>> {log}
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
            "logs/antismash/antismash/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            genefinding="none",
        shell:
            """
            antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} --logfile {log} 2>> {log}
            """

elif antismash_major_version >= 7:
    rule antismash_db_setup:
        output:
            directory(f"resources/antismash{antismash_major_version}_db"),
        conda:
            "../envs/antismash.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        shell:
            """
            download-antismash-databases --database-dir {output} &>> {log}
            antismash --version >> {log}
            antismash --database {output} --prepare-data &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            resources=f"resources/antismash{antismash_major_version}_db/",
        output:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
            json="data/interim/antismash/{version}/{strains}/{strains}.json",
            zip="data/interim/antismash/{version}/{strains}/{strains}.zip",
        conda:
            "../envs/antismash.yaml"
        threads: 4
        log:
            "logs/antismash/antismash/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            genefinding="none",
        shell:
            """
            antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} \
                --database {input.resources} \
                --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} --logfile {log} 2>> {log}
            """

rule copy_antismash:
    input:
        dir="data/interim/antismash/{version}/{strains}",
    output:
        dir=directory("data/processed/{name}/antismash/{version}/{strains}"),
    conda:
        "../envs/antismash.yaml"
    log:
        "logs/antismash/copy_antismash/copy_antismash_{version}-{strains}-{name}.log",
    shell:
        """
        base_dir=$PWD
        mkdir {output.dir}
        (cd {output.dir} && for item in $(ls $base_dir/{input.dir}); do ln -s $base_dir/{input.dir}/$item $(basename $item); done) 2>> {log}
        """
