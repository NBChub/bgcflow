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
    log: "workflow/report/logs/roary/roary-{name}.log"
    shell:
        """
        roary -p {threads} -f {output.roary_dir} -e -n -i {params.i} -cd {params.core} -r -v --mafft {input.gff} >> {log}
        """

rule eggnog_roary:
    input:
        faa = "data/interim/roary/{name}/pan_genome_reference.fa",
        eggnog_db = "resources/eggnog_db",
        dmnd = "resources/eggnog_db/bacteria.dmnd"
    output:
        eggnog_dir = directory("data/interim/roary/{name}/eggnog/")
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log: "workflow/report/logs/roary/eggnog-{name}.log"
    shell:
        """
        mkdir -p {output.eggnog_dir}
        emapper.py -i {input.faa} --decorate_gff "yes" --excel --cpu {threads} -o {wildcards.strains} --output_dir {output.eggnog_dir} --data_dir {input.eggnog_db} >> {log}
        """