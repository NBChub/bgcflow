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

rule bigslice:
    input: 
        resource = "resources/bigslice/install_note.txt",
        tmp_dir = "data/interim/bgcs/{name}/{version}/",
        taxonomy = "data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
    output:
        folder = directory("data/interim/bigslice/{name}_antismash_{version}/")
    conda:
        "../envs/bigslice.yaml"
    threads: 16
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
        tmp_dir = "data/interim/bgcs/{name}/{version}/",
        bigslice_dir = "resources/bigslice/full_run_result/"
    output:
        folder = directory("data/interim/bigslice/query/{name}_antismash_{version}/")
    conda:
        "../envs/bigslice.yaml"
    threads: 8
    log:
        "workflow/report/logs/bigslice/query_bigslice/query_bigslice_{name}-antismash-{version}.log"
    params:
        n_ranks = 10,
        query_name = "{name}"
    shell:
        """
        bigslice --query {input.tmp_dir} --n_ranks {params.n_ranks} {input.bigslice_dir} -t {threads} --query_name {params.query_name} 2>> {log}
        python workflow/bgcflow/bgcflow/data/get_bigslice_query_result.py {params.query_name} {output.folder} {input.bigslice_dir} 2>> {log}
        """
