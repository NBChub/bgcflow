rule clinker:
    input:
        gbk=lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
    output:
        txt = "data/processed/{name}/clinker/{version}/clinker.txt",
        html = "data/processed/{name}/clinker/{version}/clinker.html",
    conda:
        "../envs/clinker.yaml"
    threads: 2
    log: "workflow/report/logs/clinker/clinker-{version}-{name}.log"
    shell:
        """
        clinker {input.gbk} -j {threads} -o {output.txt} -p {output.html} 2>> {log}
        """
