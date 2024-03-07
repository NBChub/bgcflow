rule gecco:
    input:
        genbank = "data/interim/processed-genbank/{strains}.gbk"
    output:
        gecco = "data/interim/gecco/{gecco_version}/{strains}/{strains}.clusters.tsv",
        #antismash_sideload = "data/interim/gecco/{gecco_version}/{strains}/{strains}.sideload.json",
    conda: "../envs/gecco.yaml"
    log: "logs/gecco/run/{gecco_version}_{strains}.log"
    threads: 2
    shell:
        """
        gecco run --genome {input.genbank} \
            --output-dir data/interim/gecco/{wildcards.gecco_version}/{wildcards.strains} \
            --antismash-sideload \
            --jobs {threads} \
            --force-tsv 2>> {log}
        """

rule antismash_sideload_gecco:
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
        tigrfam=rules.antismash_db_setup.output.tigrfam,
        #gecco=rules.gecco.output.antismash_sideload
    output:
        folder=directory("data/interim/gecco/{gecco_version}/antismash_sideload/{version}/{strains}/"),
        gbk="data/interim/gecco/{gecco_version}/antismash_sideload/{version}/{strains}/{strains}.gbk",
        json="data/interim/gecco/{gecco_version}/antismash_sideload/{version}/{strains}/{strains}.json",
        zip="data/interim/gecco/{gecco_version}/antismash_sideload/{version}/{strains}/{strains}.zip",
    conda:
        "../envs/antismash.yaml"
    threads: 4
    log:
        "logs/gecco/{gecco_version}/antismash_sideload/{version}/antismash_{version}-{strains}.log",
    params:
        folder=directory("data/interim/antismash/{version}/{strains}/"),
        antismash_db_path=antismash_db_path,
        genefinding="none",
    shell:
        """
        # Find the latest existing JSON output for this strain
        latest_version=$(find data/interim/antismash/*/{wildcards.strains} -name "{wildcards.strains}.json" -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f 2- | cut -d '/' -f 4) 2>> {log}

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
        antismash --genefinding-tool {params.genefinding} \
            --output-dir {params.folder} \
            --database {params.antismash_db_path} \
            --cb-general \
            --cb-subclusters \
            --cb-knownclusters \
            -c {threads} $antismash_input \
            --sideload {input.gecco} \
            --logfile {log} 2>> {log}

        # Check if the run failed due to changed detection results
        if grep -q "ValueError: Detection results have changed. No results can be reused" {log}; then
            # Use genbank input instead
            echo "Previous JSON result is invalid, starting AntiSMASH from scratch..." >> {log}
            antismash --genefinding-tool {params.genefinding} \
                --output-dir {params.folder} \
                --database {params.antismash_db_path} \
                --cb-general \
                --cb-subclusters \
                --cb-knownclusters \
                -c {threads} {input.gbk} \
                --sideload {input.gecco} \
                --logfile {log} 2>> {log}
        fi
        """

rule gecco_aggregate:
    input:
        gecco = lambda wildcards: expand(
            "data/interim/gecco/{gecco_version}/{strains}/{strains}.clusters.tsv",
            name=wildcards.name,
            gecco_version=wildcards.gecco_version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        #antismash_sideload_gecco = lambda wildcards: expand(
        #    "data/interim/gecco/{gecco_version}/antismash_sideload/{version}/{strains}/{strains}.json",
        #    name=wildcards.name,
        #    gecco_version=wildcards.gecco_version,
        #    version=wildcards.version,
        #    strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        #)
    output:
        #gecco = directory("data/processed/{name}/gecco/{gecco_version}/antismash_sideload/{version}"),
        gecco = directory("data/processed/{name}/gecco/{gecco_version}"),
    conda: "../envs/gecco.yaml"
    log: "logs/gecco/aggregate/{name}_{gecco_version}.log"
    threads: 2
    shell:
        """
        set -e
        TMPDIR="data/interim/gecco/tmp/{wildcards.name}/{wildcards.gecco_version}"
        mkdir -p $TMPDIR
        INPUT_TSV="$TMPDIR/gecco_output.txt"
        echo '{input.gecco}' | tr ' ' '\n' > $INPUT_TSV
        python workflow/bgcflow/bgcflow/data/gecco_aggregate.py $INPUT_TSV {output.gecco}/gecco_clusters.csv 2>> {log}
        """
