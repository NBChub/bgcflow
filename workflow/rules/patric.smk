if PATRIC == []:
    pass
else:
    rule patric_genome_download:
        output:
            fna = "data/interim/fasta/{patric}.fna",
        shell:
            """
            cd data/interim/fasta
            wget -qN "ftp://ftp.patricbrc.org/genomes/{wildcards.patric}/{wildcards.patric}.fna"
            """