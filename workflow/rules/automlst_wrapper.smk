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

rule autoMLST_wrapper:
    input:
        gbk = expand("data/processed/genbank/{strains}.gbk", strains=STRAINS),
        reduced_core = "resources/automlst-simplified-wrapper-main/reducedcore.hmm"
    output:
        automlst_dir = directory("data/interim/automlst_wrapper/"),
        tree = "data/interim/automlst_wrapper/raxmlpart.txt.treefile"        
    conda: "../envs/automlst_wrapper.yaml"
    shell:
        """
        cp {input.gbk} {output.automlst_dir}
        python resources/automlst-simplified-wrapper-main/simplified_wrapper.py {output.automlst_dir}
        """