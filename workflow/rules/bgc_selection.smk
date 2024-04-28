def get_bgc_inputs(pep_object, antismash_version):
    """
    Given a PEP Object, get all genbank files based on the bgc_id column
    """
    antismash_path = Path(f"data/interim/antismash/{antismash_version}")
    gbk_list = []

    base_path = Path(pep_object.config["sample_table"]).parent
    input_path = Path(".")
    if "input_folder" in list(pep_object.config.keys()):
        input_path = base_path / pep_object.config["input_folder"]
        input_path = input_path.resolve()

    df = pep_object.sample_tables
    for i in df.index:
        bgc_id = df.loc[i, "bgc_id"]
        genome_id = df.loc[i, "genome_id"]
        # override with custom path
        if 'gbk_path' in df.columns and 'input_file' not in df.columns:
            df.rename(columns={'gbk_path': 'input_file'}, inplace=True)
        assert 'input_file' in df.columns
        custom_path = df.loc[i, "input_file"]

        if custom_path != None:
            gbk_path = Path(input_path) / custom_path
        elif 'input_folder' in pep_object.config.keys():
            gbk_path = Path(input_path / f"{bgc_id}.gbk")
        else:
            gbk_path = antismash_path / genome_id / f"{bgc_id}.gbk"
        gbk_list.append(gbk_path)
    return gbk_list

rule downstream_bgc_prep_selection:
    input:
        gbk=lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
        table="data/processed/{name}/tables/df_gtdb_meta.csv",
    output:
        input_list=temp("data/interim/bgcs/{name}/{version}/input_list.txt"),
        taxonomy="data/interim/bgcs/taxonomy/taxonomy_{name}_antismash_{version}.tsv",
        outdir=directory("data/interim/bgcs/{name}/{version}"),
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        dataset="data/interim/bgcs/datasets.tsv",
    log: "logs/bgcs/downstream_bgc_prep/{name}/downstream_bgc_prep-{version}.log",
    shell:
        """
        echo "Preparing BGCs for {wildcards.name} downstream analysis..." >> {log}

        echo "Step 1. Generate symlink for each regions in genomes in dataset" >> {log}
        echo {input.gbk} | tr ' ' '\n' > {output.input_list} 2>> {log}
        head -n 5 {output.input_list} >> {log}
        python workflow/bgcflow/bgcflow/data/bgc_downstream_prep_selection.py {output.input_list} {output.outdir} 2>> {log}

        echo "Step 2. Generate taxonomic information for dataset" >> {log}
        python workflow/bgcflow/bgcflow/data/bigslice_prep.py {input.table} {output.taxonomy} 2>> {log}
        # append new dataset information
        ## check if previous dataset exists
        if [[ -s {params.dataset} ]]
        then
            echo "Previous dataset detected, appending dataset information for {wildcards.name}..."
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        else
            echo "No previous dataset detected, generating dataset information for {wildcards.name}..." 2>> {log}
            echo -e '# Dataset name\tPath to folder\tPath to taxonomy\tDescription' > {params.dataset} 2>> {log}
            sed -i 'a {wildcards.name}_antismash_{wildcards.version}\t{wildcards.name}_antismash_{wildcards.version}\ttaxonomy/taxonomy_{wildcards.name}_antismash_{wildcards.version}.tsv\t{wildcards.name}' {params.dataset} 2>> {log}
        fi
        echo "Step 3. Generate mapping for visualization" >> {log}
        python workflow/bgcflow/bgcflow/data/get_bigscape_mapping.py {output.outdir} {output.bgc_mapping} 2>> {log}
        """
