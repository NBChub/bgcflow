rule ncbi_genome_download:
    output:
        "data/interim/fasta/{strains}.fna",
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s refseq -F fasta -A {strains} -o data/raw/ncbi/download bacteria --verbose
        gunzip -c data/raw/ncbi/download/refseq/bacteria/{strains}/*.fna.gz > {output}
        #rm -rf data/raw/ncbi/download/refseq/bacteria/{strains}
        """
