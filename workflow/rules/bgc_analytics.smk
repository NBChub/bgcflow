
rule antismash_summary:
    input: 
        antismash_dir = "data/interim/bgcs/{name}/{version}",
        final_copy = lambda wildcards: expand("data/processed/{name}/antismash/{version}/{strains}",
                                               name=wildcards.name, version=wildcards.version,
                                               strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)])
    output:
        df_antismash_summary = report("data/processed/{name}/tables/df_antismash_{version}_summary.csv", caption="../report/table-antismash.rst", category="{name}", subcategory="AntiSMASH Summary Table")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/bgc_analytics/antismash_summary-{version}-{name}.log"
    params:
        df_samples = lambda wildcards: PEP_PROJECTS[wildcards.name].config['sample_table'],
        fna_dir = "data/interim/fasta/",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_genome_dataset.py {params.fna_dir} {input.antismash_dir} '{params.df_samples}' {output.df_antismash_summary} 2> {log}
        """

rule write_dependency_versions:
    output:
        "workflow/report/dependency_versions.json"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/bgc_analytics/write_dependency_versions.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_dependencies.py {output} 2> {log}
        """
