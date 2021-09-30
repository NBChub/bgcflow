rule roary:
    input: 
        expand("data/interim/prokka/{strains}/{strains}.gff", strains = STRAINS),
    output:
        directory("data/processed/roary/")
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
    threads: 12
    shell:
        """
        roary -p {threads} -f {output} -e -n -i {params.i} -v --mafft {input}
        """