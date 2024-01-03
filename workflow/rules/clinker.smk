rule antismash_colourmap:
    output:
        colour_map = "data/interim/clinker/colours.txt"
    log:
        "logs/clinker/antismash_colourmap.log"
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        biosynthetic = "#810e15",
        biosynthetic_additional = "#f16d75",
        transport = "#6495ed",
        regulatory = "#2e8b57",
        other = "#808080",
        resistance = "#ffc24b"
    shell:
        """
        echo "biosynthetic,{params.biosynthetic}" > {output.colour_map}
        echo "biosynthetic-additional,{params.biosynthetic_additional}" >> {output.colour_map}
        echo "transport,{params.transport}" >> {output.colour_map}
        echo "regulatory,{params.regulatory}" >> {output.colour_map}
        echo "other,{params.other}" >> {output.colour_map}
        echo "resistance,{params.resistance}" >> {output.colour_map}
        """
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
rule clinker_gene_functions:
    input:
        gbk_dir = "data/interim/clinker/{name}/{version}",
    output:
        gene_functions = "data/interim/clinker/{name}/{version}__gene_functions.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/clinker/clinker_gene_functions_{name}_{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/features/get_antismash_gene_kind.py {input.gbk_dir} {output.gene_functions} 2>> {log}
        """
rule clinker:
    input:
        gbk_dir = "data/interim/clinker/{name}/{version}",
        colour_map = "data/interim/clinker/colours.txt",
        gene_functions = "data/interim/clinker/{name}/{version}__gene_functions.csv",
    output:
        txt = "data/processed/{name}/clinker/{version}/clinker.txt",
        html = "data/processed/{name}/clinker/{version}/clinker.html",
    conda:
        "../envs/clinker.yaml"
    threads: 4
    log: "logs/clinker/clinker-{version}-{name}.log"
    shell:
        """
        clinker {input.gbk_dir}/*.gbk -gf {input.gene_functions} -cm {input.colour_map} -j {threads} -o {output.txt} -p {output.html} 2>> {log}
        """

rule clinker_extract:
    input:
        gbk_dir = "data/interim/clinker/{name}/{version}",
        txt = "data/processed/{name}/clinker/{version}/clinker.txt",
    output:
        csv = "data/processed/{name}/clinker/{version}/clinker.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    threads: 4
    log: "logs/clinker/clinker_extract-{version}-{name}.log"
    params:
        cutoff = 0.7
    shell:
        """
        python workflow/bgcflow/bgcflow/features/clinker_extract.py {input.txt} {input.gbk_dir} {output.csv} {params.cutoff} 2>> {log}
        """
