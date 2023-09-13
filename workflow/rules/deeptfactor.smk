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
        faa="data/interim/prokka/{strains}/{strains}.faa",
        resource="resources/deeptfactor/",
    output:
        deeptfactor="data/interim/deeptfactor/{strains}/prediction_result.txt",
    conda:
        "../envs/deeptfactor.yaml"
    threads: 2
    params:
        outdir="data/interim/deeptfactor/{strains}/",
    log:
        "logs/deeptfactor/deeptfactor/deeptfactor-{strains}.log",
    shell:
        """
        workdir=$PWD
        mkdir -p data/interim/deeptfactor/{wildcards.strains} 2>> {log}
        (cd {input.resource} && python tf_running.py \
            -i $workdir/{input.faa} -o $workdir/{params.outdir} \
            -g cpu -cpu {threads}) 2>> {log}
        """


rule deeptfactor_to_json:
    input:
        deeptfactor="data/interim/deeptfactor/{strains}/prediction_result.txt",
    output:
        deeptfactor_json="data/interim/deeptfactor/{strains}_deeptfactor.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/deeptfactor/deeptfactor/deeptfactor-{strains}_to_json.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/deeptfactor_scatter.py {input.deeptfactor} {output.deeptfactor_json} 2>> {log}
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
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_deeptfactor.txt"
        echo '{input.deeptfactor_json}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/deeptfactor_gather.py $INPUT_JSON {output.df_deeptfactor} 2>> {log}
        rm $INPUT_JSON
        """
