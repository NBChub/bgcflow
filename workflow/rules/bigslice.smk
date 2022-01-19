rule install_bigslice:
    output:
        "resources/bigslice.txt"
    conda:
        "../envs/bigslice.yaml"
    shell:
        """
        (cd resources && download_bigslice_hmmdb && rm bigslice_models.tar.gz)
        bigslice --version . > {output}
        """

rule get_bigslice_inputs:
    input:
        gbk = lambda wildcards: get_bigscape_inputs(wildcards.name, wildcards.version),
        table = "data/processed/tables/df_gtdb_meta-{name}.csv"
    output:
        taxonomy = "data/interim/bigslice/tmp/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        tempdir = directory("data/interim/bigslice/tmp/{name}_antismash_{version}")
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/bigslice/{name}_antismash_{version}/temp_bigscape.log"
    params:
        dataset = "data/interim/bigslice/tmp/datasets.tsv"
    shell:
        """
        mkdir -p {output.tempdir} 2> {log}
        for i in $(dirname {input.gbk})
        do
            for r in $(ls $i/*.region*.gbk)
            do
                parent_dir=$(dirname $PWD/$r)
                filename=$(basename $r)
                mkdir {output.tempdir}/$(basename $parent_dir)
                python workflow/bgcflow/bgcflow/data/bigslice_prep.py {input.table} {output.taxonomy}
                (cd {output.tempdir}/$(basename $parent_dir) && ln -s $parent_dir/$filename $filename) 2>> {log}
            done
        echo -e '# Dataset name\tPath to folder\tPath to taxonomy\tDescription' > {params.dataset}
        sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset}
        done
        """

rule bigslice:
    input: 
        ref = "resources/bigslice.txt",
        dir = directory("data/interim/bigslice/tmp/{name}_antismash_{version}/"),
        taxonomy = "data/interim/bigslice/tmp/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
    output:
        log ="data/interim/bigslice/{name}_antismash_{version}_check.txt"
    conda:
        "../envs/bigslice.yaml"
    threads: 64
    params:
        n_ranks = 10,
        folder = "data/interim/bigslice/{name}_antismash_{version}/",
    shell:
        """
        bigslice -i data/interim/bigslice/tmp/ {params.folder}
        echo finish > {output.log}
        """
