rule mkdocs_py_report:
    output:
        notebook = "data/processed/{name}/docs/{bgcflow_rules}.ipynb",
        markdown = "data/processed/{name}/docs/{bgcflow_rules}.md",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/{bgcflow_rules}/{bgcflow_rules}-report-{name}.log"
    wildcard_constraints:
        name="|".join(py_wildcards),
    params:
        notebook = "workflow/notebook/{bgcflow_rules}.py.ipynb"
    shell:
        """
        cp {params.notebook} {output.notebook}
        jupyter nbconvert --to notebook --execute {output.notebook} --output {wildcards.bgcflow_rules}.ipynb 2>> {log}
        jupyter nbconvert --to markdown {output.notebook} --no-input --output {wildcards.bgcflow_rules}.md 2>> {log}
        """

rule mkdocs_rpy_report:
    output:
        notebook = "data/processed/{name}/docs/{bgcflow_rules}.ipynb",
        markdown = "data/processed/{name}/docs/{bgcflow_rules}.md",
    conda:
        "../envs/r_notebook.yaml"
    log: "workflow/report/logs/{bgcflow_rules}/{bgcflow_rules}-report-{name}.log"
    wildcard_constraints:
        name="|".join(rpy_wildcards),
    params:
        notebook = "workflow/notebook/{bgcflow_rules}.rpy.ipynb"
    shell:
        """
        cp {params.notebook} {output.notebook}
        jupyter nbconvert --to notebook --execute {output.notebook} --output {wildcards.bgcflow_rules}.ipynb 2>> {log}
        jupyter nbconvert --to markdown {output.notebook} --no-input --output {wildcards.bgcflow_rules}.md 2>> {log}
        """