antismash_db_path = Path(f"resources/antismash_db")

if antismash_major_version <= 6:
    rule antismash_db_setup:
        output:
            clusterblast=directory(antismash_db_path / "clusterblast"),
            clustercompare_mibig=directory(antismash_db_path / "clustercompare/mibig"),
            pfam=directory(antismash_db_path / "pfam/34.0"),
            resfam=directory(antismash_db_path / "resfam"),
            tigrfam=directory(antismash_db_path / "tigrfam"),
        conda:
            "../envs/antismash_v6.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        params:
            antismash_db_path=antismash_db_path
        shell:
            """
            download-antismash-databases --database-dir {params.antismash_db_path} &>> {log}
            antismash --version >> {log}
            antismash --check-prereqs --databases {params.antismash_db_path} &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            clusterblast=rules.antismash_db_setup.output.clusterblast,
            clustercompare_mibig=rules.antismash_db_setup.output.clustercompare_mibig,
            pfam=rules.antismash_db_setup.output.pfam,
            resfam=rules.antismash_db_setup.output.resfam,
            tigrfam=rules.antismash_db_setup.output.tigrfam
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
            antismash_db_path=antismash_db_path,
        shell:
            """
            antismash \
                --genefinding-tool {params.genefinding} \
                --database {params.antismash_db_path} \
                --output-dir {params.folder} \
                --cb-general \
                --cb-subclusters \
                --cb-knownclusters \
                -c {threads} {input.gbk} --logfile {log} 2>> {log}
            """

elif antismash_major_version >= 7:
    rule antismash_db_setup:
        output:
            as_js=directory(antismash_db_path / "as-js/0.13"),
            clusterblast=directory(antismash_db_path / "clusterblast"),
            clustercompare_mibig=directory(antismash_db_path / "clustercompare/mibig/3.1"),
            comparippson_asdb=directory(antismash_db_path / "comparippson/asdb/4.0"),
            comparippson_mibig=directory(antismash_db_path / "comparippson/mibig/3.1"),
            knownclusterblast=directory(antismash_db_path / "knownclusterblast/3.1"),
            nrps_pks_stachelhaus=directory(antismash_db_path / "nrps_pks/stachelhaus/1.1"),
            nrps_pks_svm=directory(antismash_db_path / "nrps_pks/svm/2.0"),
            nrps_pks_transATor=directory(antismash_db_path / "nrps_pks/transATor/2023.02.23"),
            pfam=directory(antismash_db_path / "pfam/35.0"),
            resfam=directory(antismash_db_path / "resfam"),
            tigrfam=directory(antismash_db_path / "tigrfam"),
        conda:
            "../envs/antismash.yaml"
        log:
            "logs/antismash/antismash_db_setup.log",
        params:
            antismash_db_path=antismash_db_path
        shell:
            """
            download-antismash-databases --database-dir {params.antismash_db_path} &>> {log}
            antismash --version >> {log}
            antismash --database {params.antismash_db_path} --prepare-data &>> {log}
            """

    rule antismash:
        input:
            gbk="data/interim/processed-genbank/{strains}.gbk",
            as_js=rules.antismash_db_setup.output.as_js,
            clusterblast=rules.antismash_db_setup.output.clusterblast,
            clustercompare_mibig=rules.antismash_db_setup.output.clustercompare_mibig,
            comparippson_asdb=rules.antismash_db_setup.output.comparippson_asdb,
            comparippson_mibig=rules.antismash_db_setup.output.comparippson_mibig,
            knownclusterblast=rules.antismash_db_setup.output.knownclusterblast,
            nrps_pks_stachelhaus=rules.antismash_db_setup.output.nrps_pks_stachelhaus,
            nrps_pks_svm=rules.antismash_db_setup.output.nrps_pks_svm,
            nrps_pks_transATor=rules.antismash_db_setup.output.nrps_pks_transATor,
            pfam=rules.antismash_db_setup.output.pfam,
            resfam=rules.antismash_db_setup.output.resfam,
            tigrfam=rules.antismash_db_setup.output.tigrfam
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
            antismash_db_path=antismash_db_path,
            genefinding="none",
            taxon="bacteria",
        shell:
            """
            set +e

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
            antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} \
                --database {params.antismash_db_path} --taxon {params.taxon} \
                --cb-general --cb-subclusters --cb-knownclusters -c {threads} $antismash_input --logfile {log} 2>> {log}

            # Check if the run failed due to changed detection results or changed protocluster types
            if grep -q -e "ValueError: Detection results have changed. No results can be reused" \
                    -e "RuntimeError: Protocluster types supported by HMM detection have changed, all results invalid" {log}
            then
                # Use genbank input instead
                echo "Previous JSON result is invalid, starting AntiSMASH from scratch..." >> {log}
                antismash --genefinding-tool {params.genefinding} --output-dir {params.folder} \
                    --database {params.antismash_db_path} --taxon {params.taxon} \
                    --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input.gbk} --logfile {log} 2>> {log}
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
