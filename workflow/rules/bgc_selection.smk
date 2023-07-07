def get_bgc_inputs(pep_object, antismash_version):
    """
    Given a PEP Object, get all genbank files based on the bgc_id column
    """
    antismash_path = Path(f"data/interim/antismash/{antismash_version}")
    gbk_list = []

    base_path = Path(pep_object.config["sample_table"]).parent
    input_path = ""
    if "input_folder" in list(pep_object.config.keys()):
        input_path = base_path / pep_object.config["input_folder"]
        input_path = input_path.resolve()

    df = pep_object.sample_tables
    for i in df.index:
        bgc_id = df.loc[i, "bgc_id"]
        genome_id = df.loc[i, "genome_id"]
        # override with custom path
        assert 'gbk_path' in df.columns
        custom_path = df.loc[i, "gbk_path"]
        #print(custom_path, type(custom_path), custom_path != None, file=sys.stderr)

        if custom_path != None:
            gbk_path = custom_path
        elif 'input_folder' in pep_object.config.keys():
            gbk_path = Path(input_path / f"{bgc_id}.gbk")
        else:
            gbk_path = antismash_path / genome_id / f"{bgc_id}.gbk"
        #print(bgc_id, gbk_path, file=sys.stderr)
        gbk_list.append(gbk_path)
    return gbk_list

rule downstream_bgc_prep_selection:
    input:
        gbk=lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
        table="data/processed/{name}/tables/df_gtdb_meta.csv",
    output:
        taxonomy="data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        outdir=directory("data/interim/bgcs/{name}/{version}"),
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        dataset="data/interim/bgcs/datasets.tsv",
    log:
        general="logs/bgcs/downstream_bgc_prep/{name}/downstream_bgc_prep-{version}.log",
        symlink="logs/bgcs/downstream_bgc_prep/{name}/bgc_downstream_bgc_prep-{version}.log",
        taxonomy="logs/bgcs/downstream_bgc_prep/{name}/tax_downstream_bgc_prep-{version}.log",
    shell:
        """
        echo "Preparing BGCs for {wildcards.name} downstream analysis..." 2>> {log.general}
        #mkdir -p {output.outdir} 2>> {log.general}
        # Generate symlink for each regions in genomes in dataset
        for gbk in {input.gbk}
        do
            echo $gbk $(dirname $gbk) 2>> {log.symlink}
            i=$(dirname $gbk)
            echo Processing $i 2>> {log.symlink}
            python workflow/bgcflow/bgcflow/data/bgc_downstream_prep_selection.py $i {output.outdir} '{input.gbk}' 2>> {log.symlink}
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
