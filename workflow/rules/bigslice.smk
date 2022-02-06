rule install_bigslice:
    output:
        "resources/bigslice/install_note.txt"
    conda:
        "../envs/bigslice.yaml"
    log:
        "workflow/report/logs/bigslice/install_bigslice.log"
    shell:
        """
        (cd resources && download_bigslice_hmmdb && rm bigslice_models.tar.gz) 2>> {log}
        bigslice --version . > {output} 2>> {log}
        """

rule get_bigslice_inputs:
    input:
        gbk = lambda wildcards: get_antismash_inputs(wildcards.name, wildcards.version),
        table = "data/processed/{name}/tables/df_gtdb_meta.csv"
    output:
        taxonomy = "data/interim/bigslice/tmp/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        tempdir = directory("data/interim/bigslice/tmp/{name}_antismash_{version}")
    conda:
        "../envs/bigslice.yaml"
    log: "workflow/report/logs/bigslice/get_bigslice_inputs/get_bigslice_inputs-{name}_antismash_{version}.log"
    params:
        dataset = "data/interim/bigslice/tmp/datasets.tsv"
    shell:
        """
        echo "Preparing BiG-SLICE input for {wildcards.name}..." 2>> {log}
        mkdir -p {output.tempdir} 2>> {log}

        # Generate symlink for each regions in genomes in dataset
        for i in $(dirname {input.gbk})
        do
            mkdir {output.tempdir}/$(basename $i) 2>> {log}
            for r in $(ls $i/*.region*.gbk)
            do
                parent_dir=$(dirname $PWD/$r)
                filename=$(basename $r)
                (cd {output.tempdir}/$(basename $parent_dir) && ln -sf $parent_dir/$filename $filename) 2>> {log}
            done
        done

        # generate taxonomic information for dataset
        python workflow/bgcflow/bgcflow/data/bigslice_prep.py {input.table} {output.taxonomy} 2>> {log}

        # append new dataset information
        ## check if previous dataset exists
        if [[ -s {params.dataset} ]]
        then
            echo "Previous dataset detected, appending dataset information for {wildcards.name}..."
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        else
            echo "No previous dataset detected, generating dataset information for {wildcards.name}..."
            echo -e '# Dataset name\tPath to folder\tPath to taxonomy\tDescription' > {params.dataset} 2>> {log}
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        fi
        """

rule bigslice:
    input: 
        resource = "resources/bigslice/install_note.txt",
        tmp_dir = "data/interim/bigslice/tmp/{name}_antismash_{version}/",
        taxonomy = "data/interim/bigslice/tmp/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
    output:
        folder = directory("data/interim/bigslice/{name}_antismash_{version}/")
    conda:
        "../envs/bigslice.yaml"
    threads: 64
    log:
        "workflow/report/logs/bigslice/bigslice/bigslice_{name}-antismash-{version}.log"
    shell:
        """
        bigslice -i data/interim/bigslice/tmp/ {output.folder} -t {threads} > {log} 2>> {log}
        """

rule fetch_bigslice_db:
    input:
        resource = "resources/bigslice/install_note.txt",
    output:
        folder = directory("resources/bigslice/full_run_result/")
    conda:
        "../envs/bigslice.yaml"
    log:
        "workflow/report/logs/bigslice/fetch_bigslice_db.log"
    shell:
        """
        (cd resources/bigslice && wget -c -nc http://bioinformatics.nl/~kauts001/ltr/bigslice/paper_data/data/full_run_result.zip && unzip full_run_result.zip) 2>> {log}
        rm resources/bigslice/full_run_result.zip
        """

rule query_bigslice:
    input:
        tmp_dir = "data/interim/bigslice/tmp/{name}_antismash_{version}/",
        bigslice_dir = "resources/bigslice/full_run_result/"
    output:
        log = "data/interim/bigslice/query/{name}_antismash_{version}.txt"
    conda:
        "../envs/bigslice.yaml"
    threads: 32
    log:
        "workflow/report/logs/bigslice/query_bigslice/query_bigslice_{name}-antismash-{version}.log"
    params:
        n_ranks = 10,
        query_name = "{name}"
    shell:
        """
        bigslice --query {input.tmp_dir} --n_ranks {params.n_ranks} {input.bigslice_dir} -t {threads} --query_name {params.query_name} 2>> {log}
        """

