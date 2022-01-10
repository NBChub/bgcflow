rule roary:
    input: 
        gff = lambda wildcards: get_roary_inputs(wildcards.name)
    output:
        roary_dir = directory("data/interim/roary/{name}/"),
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
        core = 95,
    threads: 64
    log: "workflow/report/logs/roary_{name}.log"
    shell:
        """
        roary -p {threads} -f {output.roary_dir} -e -n -i {params.i} -cd {params.core} -r -v --mafft {input.gff} >> {log}
        """