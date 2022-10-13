if len(py_wildcards) > 0:
    rule mkdocs_py_report:
        output:
            notebook = "data/processed/{name}/docs/{bgcflow_rules_py}.ipynb",
            markdown = "data/processed/{name}/docs/{bgcflow_rules_py}.md",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "workflow/report/logs/report/{bgcflow_rules_py}-report-{name}.log"
        wildcard_constraints:
            bgcflow_rules_py="|".join(py_wildcards),
        params:
            notebook = "workflow/notebook/{bgcflow_rules_py}.py.ipynb"
        shell:
            """
            cp {params.notebook} {output.notebook} 2>> {log}
            jupyter nbconvert --to notebook --execute {output.notebook} --output {wildcards.bgcflow_rules_py}.ipynb 2>> {log}
            jupyter nbconvert --to markdown {output.notebook} --no-input --output {wildcards.bgcflow_rules_py}.md 2>> {log}
            """

if len(rpy_wildcards) > 0:
    rule mkdocs_rpy_report:
        output:
            notebook = "data/processed/{name}/docs/{bgcflow_rules_rpy}.ipynb",
            markdown = "data/processed/{name}/docs/{bgcflow_rules_rpy}.md",
        conda:
            "../envs/r_notebook.yaml"
        log: "workflow/report/logs/report/{bgcflow_rules_rpy}-report-{name}.log"
        wildcard_constraints:
            bgcflow_rules_rpy="|".join(rpy_wildcards),
        params:
            notebook = "workflow/notebook/{bgcflow_rules_rpy}.rpy.ipynb"
        shell:
            """
            cp {params.notebook} {output.notebook} 2>> {log}
            jupyter nbconvert --to notebook --execute {output.notebook} --output {wildcards.bgcflow_rules_rpy}.ipynb 2>> {log}
            jupyter nbconvert --to markdown {output.notebook} --no-input --output {wildcards.bgcflow_rules_rpy}.md 2>> {log}
            """