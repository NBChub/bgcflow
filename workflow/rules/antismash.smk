if antismash_major_version <= 6:
    rule antismash_db_setup:
        output:
            database=directory("resources/antismash_db"),
        conda:
            "../envs/antismash_v6.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        shell:
            """
            download-antismash-databases --database-dir {output.database} &>> {log}
            antismash --version >> {log}
            antismash --check-prereqs --databases {output.database} &>> {log}
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
            "../envs/antismash_v6.yaml"
        threads: 4
        log:
            "logs/antismash/antismash/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            genefinding="none",
        shell:
            """
            antismash \
                --genefinding-tool {params.genefinding} \
                --database {input.resources} \
                --output-dir {params.folder} \
                --cb-general \
                --cb-subclusters \
                --cb-knownclusters \
                -c {threads} {input.gbk} --logfile {log} 2>> {log}
            """

elif antismash_major_version >= 7:
    rule antismash_db_setup:
        output:
            directory("resources/antismash{version}_db"),
        conda:
            "../envs/antismash.yaml"
        log:
            "logs/antismash/antismash_{version}_db_setup.log",
        shell:
            """
            download-antismash-databases --database-dir {output} &>> {log}
            antismash --version >> {log}
            antismash --database {output} --prepare-data &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            resources="resources/antismash{version}_db/",
        output:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
            json="data/interim/antismash/{version}/{strains}/{strains}.json",
            zip="data/interim/antismash/{version}/{strains}/{strains}.zip",
        conda:
            "../envs/antismash.yaml"
        threads: 4
        log:
            "logs/antismash/antismash/{version}/antismash_{version}-{strains}.log",
        params:
            folder=directory("data/interim/antismash/{version}/{strains}/"),
            genefinding="none",
            cb_knownclusters="--cb-knownclusters",
            cb_subclusters="--cb-subclusters",
            cc_mibig="--cc-mibig",
            clusterhmmer="--clusterhmmer",
            tigrfam="--tigrfam",
            pfam2go="--pfam2go",
            rre="--rre",
            asf="--asf",
            tfbs="--tfbs",
        shell:
            """
            set +e
            # Define antiSMASH command
            antismash_command="antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} \
                    --database {input.resources} \
                    {params.cb_knownclusters} {params.cb_subclusters} {params.cc_mibig} {params.clusterhmmer} {params.tigrfam} {params.pfam2go} {params.rre} {params.asf} {params.tfbs} -c {threads}"

            # Find the latest existing JSON output for this strain
            latest_version=$(ls -d data/interim/antismash/*/{wildcards.strains}/{wildcards.strains}.json | grep {wildcards.strains} | sort -r | head -n 1 | cut -d '/' -f 4) 2>> {log}

            if [ -n "$latest_version" ]; then
                # Use existing JSON result as starting point
                old_json="data/interim/antismash/$latest_version/{wildcards.strains}/{wildcards.strains}.json"
                echo "Using existing JSON from $old_json as starting point..." >> {log}
                antismash_input="--reuse-result $old_json"
            else
                # No existing JSON result found, use genbank input
                echo "No existing JSON result found, starting AntiSMASH from scratch..." >> {log}
                antismash_input="{input.gbk}"
            fi

            # Run AntiSMASH
            $antismash_command $antismash_input --logfile {log} 2>> {log}

            # Check if the run failed due to changed detection results
            if grep -q -e "ValueError: Detection results have changed. No results can be reused" -e "RuntimeError: Protocluster types supported by HMM detection have changed, all results invalid" {log}; then
                # Use genbank input instead
                antismash_input="{input.gbk}"
                echo "Previous JSON result is invalid, starting AntiSMASH from scratch..." >> {log}
                $antismash_command $antismash_input --logfile {log} 2>> {log}
            fi
            """

rule copy_antismash:
    input:
        dir="data/interim/antismash/{version}/{strains}",
    output:
        dir=directory("data/processed/{name}/antismash/{version}/{strains}"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/antismash/copy_antismash/copy_antismash_{version}-{strains}-{name}.log",
    shell:
        """
        base_dir=$PWD
        mkdir {output.dir}
        (cd {output.dir} && for item in $(ls $base_dir/{input.dir}); do ln -s $base_dir/{input.dir}/$item $(basename $item); done) 2>> {log}
        """
