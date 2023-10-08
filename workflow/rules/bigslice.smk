rule bigslice_prep:
    input:
        dir="data/interim/bgcs/{name}/{version}/",
        taxonomy="data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
    output:
        folder=temp(directory("data/interim/bigslice/tmp/{name}_antismash_{version}/")),
    log:
        "logs/bigslice/bigslice_prep/bigslice_{name}-antismash-{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        project_name="{name}_antismash_{version}",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/bigslice_prep_single_project.py {input.dir} {input.taxonomy} {params.project_name} {output.folder} 2>> {log}
        """


rule bigslice:
    input:
        dir="data/interim/bigslice/tmp/{name}_antismash_{version}/",
    output:
        folder=directory("data/processed/{name}/bigslice/cluster_as_{version}/"),
    conda:
        "../envs/bigslice.yaml"
    threads: 16
    params:
        threshold=900,
    resources:
        tmpdir="data/interim/tempdir",
    log:
        "logs/bigslice/bigslice/bigslice_{name}-antismash-{version}.log",
    shell:
        """
        bigslice -i {input.dir} {output.folder} --threshold {params.threshold} -t {threads} &>> {log}
        """


rule fetch_bigslice_db:
    output:
        folder=directory("resources/bigslice/full_run_result/"),
    conda:
        "../envs/bigslice.yaml"
    log:
        "logs/bigslice/fetch_bigslice_db.log",
    shell:
        """
        (cd resources/bigslice && wget -c -nc http://bioinformatics.nl/~kauts001/ltr/bigslice/paper_data/data/full_run_result.zip && unzip full_run_result.zip) &>> {log}
        rm resources/bigslice/full_run_result.zip &>> {log}
        """


rule query_bigslice:
    input:
        tmp_dir="data/interim/bgcs/{name}/{version}/",
        bigslice_dir="resources/bigslice/full_run_result/",
    output:
        folder=directory("data/interim/bigslice/query/{name}_antismash_{version}/"),
    conda:
        "../envs/bigslice.yaml"
    threads: 8
    log:
        "logs/bigslice/query_bigslice/query_bigslice_{name}-antismash-{version}.log",
    params:
        n_ranks=10,
        query_name="{name}",
        run_id=6,
    resources:
        tmpdir="data/interim/tempdir",
    shell:
        """
        TIMESTAMP=$(date --iso-8601=hours)
        bigslice --query {input.tmp_dir} --n_ranks {params.n_ranks} {input.bigslice_dir} -t {threads} --query_name {params.query_name}_$TIMESTAMP --run_id {params.run_id} &>> {log}
        python workflow/bgcflow/bgcflow/data/get_bigslice_query_result.py {params.query_name} {output.folder} {input.bigslice_dir} &>> {log}
        """


rule summarize_bigslice_query:
    input:
        query_dir="data/interim/bigslice/query/{name}_antismash_{version}/",
    output:
        gcf_summary_csv="data/processed/{name}/bigslice/query_as_{version}/gcf_summary.csv",
        gcf_summary_json="data/processed/{name}/bigslice/query_as_{version}/gcf_summary.json",
        query_network="data/processed/{name}/bigslice/query_as_{version}/query_network.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigslice/summarize_bigslice_query/summarize_bigslice_query_{name}-antismash-{version}.log",
    params:
        bigfam_db_path="resources/bigslice/full_run_result/result/data.db",
        cutoff=900,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/summarize_bigslice_query.py {input.query_dir} data/processed/{wildcards.name}/bigslice/query_as_{wildcards.version}/ {params.bigfam_db_path} {params.cutoff} 2>> {log}
        """


rule annotate_bigfam_hits:
    input:
        region_table="data/processed/{name}/tables/df_regions_antismash_{version}.csv",
        gcf_summary_csv="data/processed/{name}/bigslice/query_as_{version}/gcf_summary.csv",
    output:
        models="data/processed/{name}/bigslice/query_as_{version}/gcf_annotation.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigslice/summarize_bigslice_query/annotate_bigslice_query_{name}-antismash-{version}.log",
    params:
        bigfam_db_path="resources/bigslice/full_run_result/result/data.db",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/annotate_bigfam_model.py {input.gcf_summary_csv} {params.bigfam_db_path} {output.models} 2>> {log}
        """
