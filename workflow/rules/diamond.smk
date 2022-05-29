rule create_diamond_db:
    input:
        faa = lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES, ext="faa")
    output:
        faa_compiled = temp("data/interim/diamond/{name}.faa"),
        diamond = "data/processed/{name}/diamond/{name}.dmnd"
    conda:
        "../envs/antismash.yaml"
    threads: 8
    log: "workflow/report/logs/create_diamond_db/create_diamond_db_{name}.log"
    params:
        diamond = "data/processed/{name}/diamond/{name}"
    shell:
        """
        cat {input.faa} > {output.faa_compiled} 2>> {log}
        diamond makedb --in {output.faa_compiled} -d {params.diamond} -p {threads} 2>> {log}
        """
