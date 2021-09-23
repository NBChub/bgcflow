rule prepare_for_annotation:
    input:
        "data/raw/internal/fasta/{strains}.fasta"
    output:
        "data/interim/fasta/{strains}.fna" 
    shell:
        """
        cat {input} | sed "s/contig/{wildcards.strains}/" | sed "s/scaffold/{wildcards.strains}_scaf/" > {output}
        """

rule prokka_refseq_setup:
    output:
        "resources/Actinos_6species.gbff"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        (cd resources && unzip Actinos_6species.zip)
        """

rule prokka:
    input: 
        fna = "data/interim/fasta/{strains}.fna",
        refseq = "resources/Actinos_6species.gbff"
    output:
        gff = "data/interim/prokka/{strains}/{strains}.gff",
    conda:
        "../envs/prokka.yaml"
    params:
        genus = "Streptomyces",
        increment = 10, 
        evalue = "1e-05"
    threads: 12
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.strains} --force --proteins {input.refseq} --prefix {wildcards.strains} --genus {params.genus} --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
        """