if PATRIC == []:
    pass
else:
    rule patric_genome_download:
        output:
            fna = "data/interim/fasta/{patric}.fna",
        conda: 
            "../envs/prokka.yaml"
        shell:
            """
            cd data/interim/fasta
            wget -qN "ftp://ftp.patricbrc.org/genomes/{wildcards.patric}/{wildcards.patric}.fna"
            """