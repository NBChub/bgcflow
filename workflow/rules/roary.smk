rule roary:
    input: 
        expand("data/interim/prokka/{strains}/{strains}.gff", strains = STRAINS + NCBI_GENOMES),
    output:
        directory("data/processed/roary/test/")
    conda:
        "../envs/roary.yaml"
    params:
        i = 80,
    threads: 12
    shell:
        """
        roary -p {threads} -f {output} -e -n -i {params.i} -v --mafft {input}
        """