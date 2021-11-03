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

rule big_slice_prep:
    input:
        #antismash_dir = "data/interim/antismash/",
        antismash = expand("data/interim/antismash/{strains}/", strains = STRAINS)
    output:
        table = "data/interim/antismash/datasets.tsv",
    run:
        dataset = {"# Dataset name":["all"], 
               "Path to folder":["."], 
               "Path to taxonomy":["all_taxonomy.tsv"], 
               "Description":["all_samples"]}
        df = pd.DataFrame.from_dict(dataset)
        df.to_csv(output.table, sep="\t", index=False)

rule bigslice:
    input: 
        ref = "resources/bigslice.txt",
        antismash = expand("data/interim/antismash/{strains}/", strains = STRAINS),
        taxonomy = "data/interim/gtdb/all_taxonomy.tsv",
    output:
        folder = directory("data/interim/bigslice/all/") 
    conda:
        "../envs/bigslice.yaml"
    threads: 12
    shell:
        """
        cp {input.taxonomy} "data/interim/antismash/all_taxonomy.tsv"
        bigslice -i data/interim/antismash/ {output.folder}
        """