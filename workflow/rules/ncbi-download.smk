rule ncbi_genome_download:
    output:
        "data/interim/fasta/{ncbi}.fna"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s refseq -F fasta -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose
        gunzip -c data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*.fna.gz > {output}
        rm -rf data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}
        """
