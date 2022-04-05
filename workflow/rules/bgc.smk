rule downstream_bgc_prep:
    input:
        gbk = lambda wildcards: get_antismash_inputs(wildcards.name, wildcards.version),
        table = "data/processed/{name}/tables/df_gtdb_meta.csv"
    output:
        taxonomy = "data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        outdir = directory("data/interim/bgcs/{name}/{version}"),
        bgc_mapping = "data/interim/bgcs/{name}/{name}_antismash_{version}.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        dataset = "data/interim/bgcs/datasets.tsv"
    log:
        general = "workflow/report/logs/bgcs/downstream_bgc_prep/{name}/downstream_bgc_prep-{version}.log",
        symlink = "workflow/report/logs/bgcs/downstream_bgc_prep/{name}/bgc_downstream_bgc_prep-{version}.log",
        taxonomy = "workflow/report/logs/bgcs/downstream_bgc_prep/{name}/tax_downstream_bgc_prep-{version}.log",
    shell:
        """
        echo "Preparing BGCs for {wildcards.name} downstream analysis..." 2>> {log.general}
        #mkdir -p {output.outdir} 2>> {log.general}
        # Generate symlink for each regions in genomes in dataset
        for i in $(dirname {input.gbk})
        do
            echo $i
            python workflow/bgcflow/bgcflow/data/bgc_downstream_prep.py $i {output.outdir} 2>> {log.symlink}
        done
        # generate taxonomic information for dataset
        python workflow/bgcflow/bgcflow/data/bigslice_prep.py {input.table} {output.taxonomy} 2>> {log.taxonomy}
        # append new dataset information
        ## check if previous dataset exists
        if [[ -s {params.dataset} ]]
        then
            echo "Previous dataset detected, appending dataset information for {wildcards.name}..."
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log.general}
        else
            echo "No previous dataset detected, generating dataset information for {wildcards.name}..."
            echo -e '# Dataset name\tPath to folder\tPath to taxonomy\tDescription' > {params.dataset} 2>> {log.general}
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log.general}
        fi
        # generate mapping for visualization
        python workflow/bgcflow/bgcflow/data/get_bigscape_mapping.py {output.outdir} {output.bgc_mapping} 2>> {log.general}
        """