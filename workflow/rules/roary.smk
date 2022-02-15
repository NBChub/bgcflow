rule roary:
    input: 
        gff = lambda wildcards: get_roary_inputs(wildcards.name)
    output:
        roary_dir = directory("data/interim/roary/{name}/"),
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
    threads: 8
    log: "workflow/report/logs/roary/roary-{name}.log"
    shell:
        """
        roary -p {threads} -f {output.roary_dir} -e -n -i {params.i} -r -v --mafft -z {input.gff} >> {log}
        """