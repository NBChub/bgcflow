rule getphylo_prep:
    input:
        gbk=lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
    output:
        getphylo_input_dir = temp(directory("data/interim/getphylo/tmp/{name}_{version}"))
    log:
        "logs/getphylo/getphylo_prep_{name}_{version}.log",
    conda:
        "../envs/getphylo.yaml"
    shell:
        """
        mkdir -p {output.getphylo_input_dir} 2>> {log}
        for file in {input.gbk}; do
            cp "$file" {output.getphylo_input_dir}/ 2>> {log}
        done

        echo "All files copied successfully!" >> {log}
        """


rule getphylo:
    input:
        gbk="data/interim/getphylo/tmp/{name}_{version}",
    output:
        getphylo_dir = directory("data/processed/{name}/getphylo/as_{version}")
    log:
        "logs/getphylo/getphylo_{name}_{version}.log",
    conda:
        "../envs/getphylo.yaml"
    threads: 8
    shell:
        """
        getphylo -g '{input.gbk}/*.gbk' -o {output.getphylo_dir} -c {threads} 2>> {log}
        """
