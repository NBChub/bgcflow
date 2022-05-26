rule arts:
    input: 
        antismash = "data/interim/antismash/{version}/{strains}/{strains}.gbk",
        resources = "resources/arts/"
    output:
        folder = directory("data/interim/arts/antismash-{version}/{strains}/"),
    conda:
        "../envs/arts.yaml"
    threads: 8
    log: "workflow/report/logs/arts/arts/arts-{version}-{strains}.log"
    params:
        ref = "resources/arts/reference/actinobacteria/"
    shell:
        """
        mkdir -p data/interim/arts/antismash-{wildcards.version}/{wildcards.strains}
        python {input.resources}/artspipeline1.py {input.antismash} {params.ref} -rd data/interim/arts/antismash-{wildcards.version}/{wildcards.strains} -cpu {threads} -opt kres,phyl 2>> {log}
        """