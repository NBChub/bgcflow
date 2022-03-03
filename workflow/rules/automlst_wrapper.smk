rule install_automlst_wrapper:
    output:
        directory("resources/automlst-simplified-wrapper-main"),
        reduced_core = "resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    conda:
        "../envs/automlst_wrapper.yaml"
    log: "workflow/report/logs/automlst_wrapper/install_automlst_wrapper.log"
    shell:
        """
        (cd resources && wget https://github.com/KatSteinke/automlst-simplified-wrapper/archive/main.zip)
        (cd resources && unzip main.zip && rm main.zip)
        (cd resources/automlst-simplified-wrapper-main && unzip reducedcore.zip) 
        """       

rule prep_automlst_gbk:
    input:
        gbk = "data/interim/processed-genbank/{strains}.gbk"
    output:
        auto_gbk = "data/interim/automlst_wrapper/{name}/{strains}.gbk",
    conda:
        "../envs/automlst_wrapper.yaml"
    log: "workflow/report/logs/automlst_wrapper/prep_automlst_gbk/prep_automlst_gbk-{name}_{strains}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/features/prep_automlst.py {input.gbk} {output.auto_gbk} {wildcards.strains}
        """

rule automlst_wrapper:
    input:
        gbk = lambda wildcards: get_automlst_inputs(wildcards.name),
        reduced_core = "resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    output:
        tree = "data/interim/automlst_wrapper/{name}/raxmlpart.txt.treefile",
    log: "workflow/report/logs/automlst_wrapper/automlst_wrapper/automlst_wrapper-{name}.log"
    conda: 
        "../envs/automlst_wrapper.yaml"
    shell:
        """
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py data/interim/automlst_wrapper/{wildcards.name}
        """

rule automlst_wrapper_out:
    input:
        tree = "data/interim/automlst_wrapper/{name}/raxmlpart.txt.treefile",
        automlst_interim = "data/interim/automlst_wrapper/{name}/",
        prokka_interim = "data/interim/prokka",
        gtdb_interim = "data/interim/gtdb",
    output:
        automlst_processed = directory("data/processed/{name}/automlst_wrapper/"),
        final_tree = "data/processed/{name}/automlst_wrapper/final.newick",
    log: "workflow/report/logs/automlst_wrapper/automlst_wrapper/automlst_wrapper_out-{name}.log"
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_phylo_tree.py {input.automlst_interim} {output.automlst_processed} {input.prokka_interim} {input.gtdb_interim} 2> {log}
        """
