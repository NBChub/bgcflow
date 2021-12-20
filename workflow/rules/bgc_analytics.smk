rule antismash_summary:
    input: 
        _all_ = expand("data/interim/antismash/{version}/{strains}/{strains}.gbk", strains=STRAINS, version=dependency_version["antismash"]),
        antismash_dir = "data/interim/antismash/",
        fna_dir = "data/interim/fasta/",
        df_samples = config["samples"],
    output:
        df_antismash_summary = "data/processed/tables/df_antismash_{version}_summary.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/antismash_{version}_summary.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_genome_dataset.py {input.fna_dir} {input.antismash_dir}/{wildcards.version} {input.df_samples} {output.df_antismash_summary} 2> {log}
        """