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
        final_tree = "data/processed/{name}/automlst_wrapper/{name}.newick"
    log: "workflow/report/logs/automlst_wrapper/automlst_wrapper/automlst_wrapper-{name}.log"
    conda: 
        "../envs/automlst_wrapper.yaml"
    shell:
        """
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py data/interim/automlst_wrapper/{wildcards.name}
        cp {output.tree} {output.final_tree}
        """
