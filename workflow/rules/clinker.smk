rule clinker:
    input: 
        antismash_dir = "data/interim/bgcs/{name}/{version}",
    output:
        txt = "data/interim/clinker/{version}/{name}/{name}_{version}.txt",
        html = "data/interim/clinker/{version}/{name}/{name}_{version}.html",
    conda:
        "../envs/clinker.yaml"
    threads: 2
    log: "workflow/report/logs/clinker/clinker-{version}-{name}.log"
    shell:
        """
        clinker {input.antismash_dir}/*/*.gbk -j {threads} -o {output.txt} -p {output.html} 2>> {log}
        """