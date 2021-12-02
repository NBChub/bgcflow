rule install_automlst_wrapper:
    output:
        directory("resources/automlst-simplified-wrapper-main"),
        reduced_core = "resources/automlst-simplified-wrapper-main/reducedcore.hmm",
    conda:
        "../envs/automlst_wrapper.yaml"
    shell:
        """
        (cd resources && wget https://github.com/KatSteinke/automlst-simplified-wrapper/archive/main.zip)
        (cd resources && unzip main.zip && rm main.zip)
        (cd resources/automlst-simplified-wrapper-main && unzip reducedcore.zip) 
        """       

rule prep_automlst_gbk:
    input:
        gbk = "data/processed/genbank/{strains}.gbk",
    output:
        auto_gbk = "data/interim/automlst_wrapper/{strains}.gbk",
    conda:
        "../envs/automlst_wrapper.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/features/prep_automlst.py {input.gbk} {output.auto_gbk} {wildcards.strains}
        """

rule automlst_wrapper:
    input:
        gbk = expand("data/interim/automlst_wrapper/{strains}.gbk", strains=STRAINS),
        reduced_core = "resources/automlst-simplified-wrapper-main/reducedcore.hmm",
        automlst_dir = "data/interim/automlst_wrapper/"
    output:
        tree = "data/interim/automlst_wrapper/raxmlpart.txt.treefile"        
    conda: 
        "../envs/automlst_wrapper.yaml"
    shell:
        """
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py {input.automlst_dir}
        """
