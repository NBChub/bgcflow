rule roary:
    input: 
        expand("data/interim/prokka/{strains}/{strains}.gff", strains = STRAINS),
    output:
        directory("data/interim/roary/")
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
    threads: 12
    shell:
        """
        roary -p {threads} -f data/interim/roary/ -e -n -i {params.i} -v --mafft {input}
        """