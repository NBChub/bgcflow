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

rule roary_out:
    input: 
        roary_interim_dir = "data/interim/roary/{name}/",
    output:
        roary_processed_dir = directory("data/processed/{name}/roary"),
        gene_presence = directory("data/processed/{name}/roary/gene_presence_absence.csv"),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/roary/roary-out-{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_pangenome_dataset.py {input.roary_interim_dir} {output.roary_processed_dir} 2> {log}
        """
