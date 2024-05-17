rule downstream_bgc_prep:
    input:
        gbk=lambda wildcards: expand("data/interim/antismash/{version}/{strains}/{strains}.gbk",
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        table="data/processed/{name}/tables/df_gtdb_meta.csv",
    output:
        input_list="data/interim/bgcs/{name}/{version}/input_list.txt",
        taxonomy="data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        outdir=directory("data/interim/bgcs/{name}/{version}"),
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        dataset="data/interim/bgcs/datasets.tsv",
        regions=lambda wildcards: get_antismash_regions(
            wildcards.name, wildcards.version, DF_SAMPLES
        ),
    log: "logs/bgcs/downstream_bgc_prep/{name}/downstream_bgc_prep-{version}.log",
    shell:
        """
        echo "Preparing BGCs for {wildcards.name} downstream analysis..." >> {log}

        echo "Step 1. Generate symlink for each regions in genomes in dataset" >> {log}
        echo {params.regions} | tr ' ' '\n' >> {output.input_list} 2>> {log}
        python workflow/bgcflow/bgcflow/data/bgc_downstream_prep_selection.py {output.input_list} {output.outdir} 2>> {log}

        echo "Step 2. Generate taxonomic information for dataset" >> {log}
        python workflow/bgcflow/bgcflow/data/bigslice_prep.py {input.table} {output.taxonomy} 2>> {log}
        # append new dataset information
        ## check if previous dataset exists
        if [[ -s {params.dataset} ]]
        then
            echo "Previous dataset detected, appending dataset information for {wildcards.name}..." >> {log}
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        else
            echo "No previous dataset detected, generating dataset information for {wildcards.name}..." >> {log}
            echo -e '# Dataset name\tPath to folder\tPath to taxonomy\tDescription' > {params.dataset} 2>> {log}
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        fi
        echo "Step 3. Generate mapping for visualization" >> {log}
        python workflow/bgcflow/bgcflow/data/get_bigscape_mapping.py {output.outdir} {output.bgc_mapping} 2>> {log}
        """
