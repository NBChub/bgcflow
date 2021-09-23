rule roary:
    input: 
        gff = expand("data/interim/prokka/{strains}/{strains}.gff", strains = STRAINS)
    output:
        "data/processed/roary/test/"
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
        tag = "matmal_dyrehaven"
    threads: 12
    shell:
        """
        roary -p {threads} -o {params.tag} -f {output} -e -n -i {params.i} -v --maft {input}
        """