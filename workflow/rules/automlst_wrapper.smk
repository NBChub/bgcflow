rule install_automlst_wrapper:
    output:
        folder=directory("resources/automlst-simplified-wrapper-main"),
        reduced_core="resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    conda:
        "../envs/automlst_wrapper.yaml"
    log:
        "workflow/report/logs/automlst_wrapper/install_automlst_wrapper.log",
    params:
        source="https://github.com/matinnuhamunada/automlst-simplified-wrapper",
        version="0.1.1"
    shell:
        """
        set -e
        mkdir -p resources 2>> {log}
        wget {params.source}/archive/refs/tags/v{params.version}.zip -O resources/automlst-simplified-wrapper-v{params.version}.zip 2>> {log}
        (cd resources && unzip -o automlst-simplified-wrapper-v{params.version}.zip && rm automlst-simplified-wrapper-v{params.version}.zip) &>> {log}
        (cd resources/automlst-simplified-wrapper-{params.version} && unzip -o reducedcore.zip) &>> {log}
        cp -r resources/automlst-simplified-wrapper-{params.version}/* {output.folder}/. 2>> {log}
        rm -rf resources/automlst-simplified-wrapper-{params.version} 2>> {log}
        """


rule prep_automlst_gbk:
    input:
        gbk="data/interim/processed-genbank/{strains}.gbk",
    output:
        auto_gbk=temp("data/interim/automlst_wrapper/{name}/{strains}.gbk"),
    conda:
        "../envs/automlst_wrapper.yaml"
    log:
        "workflow/report/logs/automlst_wrapper/prep_automlst_gbk/prep_automlst_gbk-{name}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/features/prep_automlst.py {input.gbk} {output.auto_gbk} {wildcards.strains} 2>> {log}
        """


rule automlst_wrapper:
    input:
        gbk=lambda wildcards: get_automlst_inputs(wildcards.name, DF_SAMPLES),
        reduced_core="resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    output:
        tree="data/interim/automlst_wrapper/{name}/raxmlpart.txt.treefile",
    log:
        "workflow/report/logs/automlst_wrapper/automlst_wrapper/automlst_wrapper-{name}.log",
    conda:
        "../envs/automlst_wrapper.yaml"
    threads: 8
    resources:
        tmpdir="data/interim/automlst_wrapper/tmpdir/",
    shell:
        """
        mkdir -p "data/interim/automlst_wrapper/{wildcards.name}/singles" 2>> {log}
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py data/interim/automlst_wrapper/{wildcards.name} {threads} &>> {log}
        """


rule automlst_wrapper_out:
    input:
        tree="data/interim/automlst_wrapper/{name}/raxmlpart.txt.treefile",
        organism_info=lambda wildcards: expand("data/interim/prokka/{strains}/organism_info.txt",
                    strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        automlst_processed=directory("data/processed/{name}/automlst_wrapper/"),
        final_tree="data/processed/{name}/automlst_wrapper/final.newick",
    log:
        "workflow/report/logs/automlst_wrapper/automlst_wrapper/automlst_wrapper_out-{name}.log",
    params:
        automlst_interim=lambda wildcards: f"data/interim/automlst_wrapper/{wildcards.name}/",
        prokka_interim="data/interim/prokka",
        gtdb_interim="data/interim/gtdb",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_phylo_tree.py {params.automlst_interim} {output.automlst_processed} {params.prokka_interim} {params.gtdb_interim} 2>> {log}
        """
