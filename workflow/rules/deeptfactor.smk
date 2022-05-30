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
        resource = "resources/deeptfactor/",
    output:
        deeptfactor_dir = directory("data/interim/deeptfactor/{strains}/"),
    conda:
        "../envs/deeptfactor.yaml"
    threads: 8
    params:
        faa = "../../data/interim/prokka/{strains}/{strains}.faa",
        outdir = "../../data/interim/deeptfactor/{strains}/",
    log: "workflow/report/logs/deeptfactor/deeptfactor/deeptfactor-{strains}.log"
    shell:
        """
        mkdir -p data/interim/deeptfactor/{wildcards.strains} 2>> {log}
        (cd {input.resource} && python tf_running.py \
            -i {params.faa} -o {params.outdir} \
            -g cpu -cpu {threads}) 2>> {log}
        """
