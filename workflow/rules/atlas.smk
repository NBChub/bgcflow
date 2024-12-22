rule query_atlas:
    input:
        tmp_dir="data/interim/bgcs/{name}/{version}/",
        bigslice_dir="resources/atlas",
    output:
        folder=directory("data/interim/atlas/query/{name}_antismash_{version}/"),
    conda:
        "../envs/bigslice.yaml"
    threads: 8
    log:
        "logs/bigslice/query_atlas/query_atlas_{name}-antismash-{version}.log",
    params:
        n_ranks=10,
        query_name="{name}",
        run_id=1,
        threshold=0.4,
        normalize="",
    resources:
        tmpdir="data/interim/tempdir",
    shell:
        """
        TIMESTAMP=$(date --iso-8601=hours)
        bigslice --query {input.tmp_dir} --n_ranks {params.n_ranks} {input.bigslice_dir} -t {threads} --query_name {params.query_name}_$TIMESTAMP --run_id {params.run_id} --threshold {params.threshold} {params.normalize} &>> {log}
        python workflow/bgcflow/bgcflow/data/get_bigslice_query_result.py {params.query_name} {output.folder} {input.bigslice_dir} &>> {log}
        """

rule summarize_atlas_query:
    input:
        query_dir=rules.query_atlas.output.folder,
    output:
        gcf_summary_csv="data/processed/{name}/atlas/query_as_{version}/gcf_summary.csv",
        gcf_summary_json="data/processed/{name}/atlas/query_as_{version}/gcf_summary.json",
        query_network="data/processed/{name}/atlas/query_as_{version}/query_network.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigslice/summarize_atlas_query/summarize_bigslice_query_{name}-antismash-{version}.log",
    params:
        atlas_db_path="resources/atlas/result/data.db",
        cutoff=0.4,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/summarize_bigslice_query.py {input.query_dir} data/processed/{wildcards.name}/atlas/query_as_{wildcards.version}/ {params.atlas_db_path} {params.cutoff} 2>> {log}
        """

rule annotate_atlas_hits:
    input:
        gcf_summary_csv=rules.summarize_atlas_query.output.gcf_summary_csv,
    output:
        models="data/processed/{name}/atlas/query_as_{version}/gcf_annotation.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigslice/summarize_atlas_query/annotate_bigslice_query_{name}-antismash-{version}.log",
    params:
        atlas_db_path="resources/atlas/result/data.db",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/annotate_bigfam_model.py {input.gcf_summary_csv} {params.atlas_db_path} {output.models} 2>> {log}
        """
