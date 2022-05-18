rule roary:
    input: 
        gff = lambda wildcards: get_prokka_outputs(wildcards.name)
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
        roary -p {threads} -f {output.roary_dir} -i {params.i} -e -n -r -v {input.gff} 2>> {log}
        """

rule eggnog_roary:
    input:
        roary_dir = "data/interim/roary/{name}/",
        eggnog_db = "resources/eggnog_db",
        dmnd = "resources/eggnog_db/bacteria.dmnd"
    output:
        eggnog_dir = directory("data/interim/eggnog_roary/{name}/")
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log: "workflow/report/logs/eggnog-roary/eggnog-{name}.log"
    params:
        faa = "data/interim/roary/{name}/pan_genome_reference.fa",
    shell:
        """
        mkdir -p {output.eggnog_dir}
        emapper.py -i {params.faa} --translate --itype "CDS" --excel --cpu {threads} -o {wildcards.name} --output_dir {output.eggnog_dir} --data_dir {input.eggnog_db} 2>> {log}
        """ 

rule roary_out:
    input: 
        roary_interim_dir = "data/interim/roary/{name}/",
        automlst_processed_dir = "data/processed/{name}/automlst_wrapper/",
    output:
        roary_processed_dir = directory("data/processed/{name}/roary"),
        gene_presence = "data/processed/{name}/roary/df_gene_presence_binary.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/roary/roary-out-{name}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_pangenome_dataset.py {input.roary_interim_dir} {output.roary_processed_dir} {input.automlst_processed_dir} 2>> {log}
        """
