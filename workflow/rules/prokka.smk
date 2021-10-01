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
    threads: 12
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.custom} --force --proteins {input.refseq} --prefix {wildcards.custom} --genus {params.genus} --strain {wildcards.custom} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
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