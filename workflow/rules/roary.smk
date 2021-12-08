rule roary:
    input: 
        gff = expand("data/interim/prokka/{strains}/{strains}.gff", strains = STRAINS),
    output:
        roary_dir = directory("data/interim/roary/all/"),
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
        core = 95,
    threads: 8
    shell:
        """
        roary -p {threads} -f {output.roary_dir} -e -n -i {params.i} -cd {params.core} -r -v --mafft {input.gff}
        """
