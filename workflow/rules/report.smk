rule mkdocs_report:
    output:
        notebook = "data/processed/{name}/docs/{bgcflow_rules}.ipynb",
        markdown = "data/processed/{name}/docs/{bgcflow_rules}.md",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/{bgcflow_rules}/{bgcflow_rules}-report-{name}.log"
    params:
        notebook = "workflow/notebook/{bgcflow_rules}.py.ipynb"
    shell:
        """
        cp {params.notebook} {output.notebook}
        jupyter nbconvert --to notebook --execute {output.notebook} --output {wildcards.bgcflow_rules}.ipynb 2>> {log}
        jupyter nbconvert --to markdown {output.notebook} --no-input --output {wildcards.bgcflow_rules}.md 2>> {log}
        """
