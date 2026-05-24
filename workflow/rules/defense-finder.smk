rule defense_finder:
    input:
        faa="data/interim/prokka/{strains}/{strains}.faa",
    output:
        defense_finder=directory("data/interim/defense_finder/{strains}/"),
    conda:
        "../envs/defense-finder.yaml"
    threads: 4
    log:
        "logs/defense-finder/defense-finder/defense-finder-{strains}.log",
    shell:
        """
        defense-finder run -w {threads} {input.faa} -o {output.defense_finder} &>> {log}
        """
