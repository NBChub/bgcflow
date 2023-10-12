rule copy_readme:
    output:
        "data/processed/{name}/README.md",
    log:
        "logs/report/copy-readme-{name}.log",

    shell:
        """
        cp workflow/notebook/README_template.md {output}
        """
if len(py_wildcards) > 0:

    rule copy_template_notebook:
        output:
            notebook="data/processed/{name}/docs/{bgcflow_rules_py}.ipynb",
        conda:
            "../envs/bgcflow_notes.yaml"
        log:
            "logs/report/{bgcflow_rules_py}-report-copy-{name}.log",
        wildcard_constraints:
            bgcflow_rules_py="|".join(py_wildcards),
        params:
            notebook="workflow/notebook/{bgcflow_rules_py}.py.ipynb",
        shell:
            """
            cp {params.notebook} {output.notebook} 2>> {log}
            """

    rule mkdocs_py_report:
        input:
            dependency_version = "data/processed/{name}/metadata/dependency_versions.json",
            notebook="data/processed/{name}/docs/{bgcflow_rules_py}.ipynb",
        output:
            markdown="data/processed/{name}/docs/{bgcflow_rules_py}.md",
        conda:
            "../envs/bgcflow_notes.yaml"
        log:
            "logs/report/{bgcflow_rules_py}-report-{name}.log",
        wildcard_constraints:
            bgcflow_rules_py="|".join(py_wildcards),
        shell:
            """
            jupyter nbconvert --to markdown --execute {input.notebook} --no-input --output {wildcards.bgcflow_rules_py}.md 2>> {log}
            """


if len(rpy_wildcards) > 0:

    rule copy_template_rnotebook:
        output:
            notebook="data/processed/{name}/docs/{bgcflow_rules_rpy}.ipynb",
        conda:
            "../envs/bgcflow_notes.yaml"
        log:
            "logs/report/{bgcflow_rules_rpy}-report-copy-{name}.log",
        wildcard_constraints:
            bgcflow_rules_rpy="|".join(rpy_wildcards),
        params:
            notebook="workflow/notebook/{bgcflow_rules_rpy}.rpy.ipynb",
        shell:
            """
            cp {params.notebook} {output.notebook} 2>> {log}
            """

    rule mkdocs_rpy_report:
        input:
            dependency_version = "data/processed/{name}/metadata/dependency_versions.json",
            notebook="data/processed/{name}/docs/{bgcflow_rules_rpy}.ipynb",
        output:
            markdown="data/processed/{name}/docs/{bgcflow_rules_rpy}.md",
        conda:
            "../envs/r_notebook.yaml"
        log:
            "logs/report/{bgcflow_rules_rpy}-report-{name}.log",
        wildcard_constraints:
            bgcflow_rules_rpy="|".join(rpy_wildcards),
        shell:
            """
            jupyter nbconvert --to markdown --execute {input.notebook} --no-input --output {wildcards.bgcflow_rules_rpy}.md 2>> {log}
            """
