rule deeptfactor_setup:
    output:
        folder = directory("resources/deeptfactor/"),
        pyfile = "resources/deeptfactor/tf_running.py",
    conda:
        "../envs/deeptfactor.yaml"
    log: "workflow/report/logs/deeptfactor/deeptfactor_setup.log"
    shell:  
        """
        git clone https://bitbucket.org/kaistsystemsbiology/deeptfactor.git {output.folder} 2>> {log}
        """

rule deeptfactor:
    input: 
        fasta = "data/interim/prokka/{strains}/{strains}.faa",
        folder = "resources/deeptfactor/",
    output:
        folder = directory("data/interim/deeptfactor/{strains}/"),
    conda:
        "../envs/deeptfactor.yaml"
    threads: 8
    log: "workflow/report/logs/deeptfactor/deeptfactor/deeptfactor-{strains}.log"
    shell:
        """
        mkdir -p data/interim/deeptfactor/{wildcards.strains} 2>> {log}
        (cd {input.folder} && python tf_running.py \
            -i "../../data/interim/prokka/{wildcards.strains}/{wildcards.strains}.faa" \
            -o "../../data/interim/deeptfactor/{wildcards.strains}/" \
            -g cpu -cpu {threads}) 2>> {log}
        """
