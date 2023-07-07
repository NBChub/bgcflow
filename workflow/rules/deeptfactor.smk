rule deeptfactor_setup:
    output:
        folder=directory("resources/deeptfactor/"),
        pyfile="resources/deeptfactor/tf_running.py",
    conda:
        "../envs/deeptfactor.yaml"
    log:
        "logs/deeptfactor/deeptfactor_setup.log",
    shell:
        """
        git clone https://bitbucket.org/kaistsystemsbiology/deeptfactor.git {output.folder} 2>> {log}
        """


rule deeptfactor:
    input:
        fasta="data/interim/prokka/{strains}/{strains}.faa",
        resource="resources/deeptfactor/",
    output:
        deeptfactor_dir=directory("data/interim/deeptfactor/{strains}/"),
    conda:
        "../envs/deeptfactor.yaml"
    threads: 2
    params:
        faa="../../data/interim/prokka/{strains}/{strains}.faa",
        outdir="../../data/interim/deeptfactor/{strains}/",
    log:
        "logs/deeptfactor/deeptfactor/deeptfactor-{strains}.log",
    shell:
        """
        mkdir -p data/interim/deeptfactor/{wildcards.strains} 2>> {log}
        (cd {input.resource} && python tf_running.py \
            -i {params.faa} -o {params.outdir} \
            -g cpu -cpu {threads}) 2>> {log}
        """


rule deeptfactor_to_json:
    input:
        deeptfactor_dir="data/interim/deeptfactor/{strains}/",
    output:
        deeptfactor_json="data/interim/deeptfactor/{strains}_deeptfactor.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/deeptfactor/deeptfactor/deeptfactor-{strains}_to_json.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/deeptfactor_scatter.py {input.deeptfactor_dir}/prediction_result.txt {output.deeptfactor_json} 2>> {log}
        """


rule deeptfactor_summary:
    input:
        deeptfactor_json=lambda wildcards: expand(
            "data/interim/deeptfactor/{strains}_deeptfactor.json",
            strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)],
        ),
    output:
        df_deeptfactor="data/processed/{name}/tables/df_deeptfactor.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/deeptfactor/deeptfactor_{name}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/deeptfactor_gather.py '{input.deeptfactor_json}' {output.df_deeptfactor} 2>> {log}
        """
