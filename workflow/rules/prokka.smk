rule prepare_for_annotation:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    shell:
        """
        cat {input} | sed "s/contig/{wildcards.custom}/" | sed "s/scaffold/{wildcards.custom}_scaf/" > {output}
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
        fna = "data/interim/fasta/{custom}.fna",
        refseq = "resources/Actinos_6species.gbff"
    output:
        gff = "data/interim/prokka/{custom}/{custom}.gff",
        gbk = "data/interim/prokka/{custom}/{custom}.gbk",
    conda:
        "../envs/prokka.yaml"
    params:
        genus = "Streptomyces",
        increment = 10, 
        evalue = "1e-05"
    threads: 8
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.custom} --force --proteins {input.refseq} --prefix {wildcards.custom} --genus {params.genus} --strain {wildcards.custom} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
        """