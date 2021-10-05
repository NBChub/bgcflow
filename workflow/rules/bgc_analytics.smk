rule antismash_summary:
    input: 
        _all_ = expand("data/interim/antismash/{strains}", strains=STRAINS),
        antismash_dir = "data/interim/antismash/",
    output:
        df_genomes = "data/processed/tables/df_genomes.csv",
        df_bgc_products = "data/processed/tables/df_bgc_products.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    script:
        "../src/data/make_genome_dataset.py"