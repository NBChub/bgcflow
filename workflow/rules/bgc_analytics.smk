rule antismash_summary:
    input: 
        _all_ = expand("data/interim/antismash/{version}/{strains}/{strains}.gbk", strains=STRAINS, version=dependency_version["antismash"]),
        antismash_dir = "data/interim/antismash/",
        fna_dir = "data/interim/fasta/",
    output:
        df_antismash_summary = report("data/processed/tables/df_antismash_{version}_summary.csv", caption="../report/table-antismash.rst", category="BGC Prediction")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/antismash_{version}_summary.log"
    params:
        df_samples = SAMPLE_PATHS,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_genome_dataset.py {input.fna_dir} {input.antismash_dir}/{wildcards.version} '{params.df_samples}' {output.df_antismash_summary} 2> {log}
        """

rule write_dependency_versions:
    output:
        "workflow/report/dependency_versions.json"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/dependency_versions.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_dependencies.py {output} 2> {log}
        """
