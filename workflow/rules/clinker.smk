rule prep_clinker:
    input:
        gbk=lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
    output:
        gbk_dir = temp(directory("data/interim/clinker/{name}/{version}"))
    log:
        "logs/clinker/prep_clinker_{name}_{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/features/prep_clinker.py '{input.gbk}' {output.gbk_dir} 2>> {log}
        """
rule clinker:
    input:
        gbk_dir = "data/interim/clinker/{name}/{version}"
    output:
        txt = "data/processed/{name}/clinker/{version}/clinker.txt",
        html = "data/processed/{name}/clinker/{version}/clinker.html",
    conda:
        "../envs/clinker.yaml"
    threads: 4
    log: "logs/clinker/clinker-{version}-{name}.log"
    shell:
        """
        clinker {input.gbk_dir}/*.gbk -j {threads} -o {output.txt} -p {output.html} 2>> {log}
        """
