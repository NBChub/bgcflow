rule gecco:
    input:
        genbank = "data/interim/processed-genbank/{strains}.gbk"
    output:
        gecco = "data/interim/gecco/{gecco_version}/{strains}/{strains}.clusters.tsv"
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

rule gecco_aggregate:
    input:
        gecco = lambda wildcards: expand(
            "data/interim/gecco/{gecco_version}/{strains}/{strains}.clusters.tsv",
            name=wildcards.name,
            gecco_version=dependency_version["gecco"],
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
    output:
        gecco = directory("data/processed/{name}/gecco/{gecco_version}")
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
        mkdir -p {output.gecco}
        while IFS= read -r file
        do
            cp "$file" {output.gecco}
        done < $INPUT_TSV
        """
