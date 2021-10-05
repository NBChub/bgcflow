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

rule prokka_ncbi:
    input: 
        fna = "data/interim/fasta/{ncbi}.fna",
        refseq = "resources/Actinos_6species.gbff"
    output:
        gff = "data/interim/prokka/{ncbi}/{ncbi}.gff",
        gbk = "data/interim/prokka/{ncbi}/{ncbi}.gbk",
        genus = temp("data/interim/prokka/{ncbi}/genus"),
        species = temp("data/interim/prokka/{ncbi}/species"),
    conda:
        "../envs/prokka.yaml"
    params:
        increment = 10, 
        evalue = "1e-05"
    threads: 12
    shell:
        """
        head -1 {input.fna} | cut -d' ' -f2 > {output.genus}
        head -1 {input.fna} | cut -d' ' -f3 > {output.species}
        if [ `cat {output.genus}` == "Streptomyces" ]
        then
            prokka --outdir data/interim/prokka/{wildcards.ncbi} --force --proteins {input.refseq} --prefix {wildcards.ncbi} --genus `cat {output.genus}` --strain {wildcards.ncbi} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
        else
            prokka --outdir data/interim/prokka/{wildcards.ncbi} --force --prefix {wildcards.ncbi} --genus `cat {output.genus}` --species `cat {output.species}` --strain {wildcards.ncbi} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input}
        fi
        """