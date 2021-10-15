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
        fna = "data/interim/fasta/{strains}.fna",
        refseq = "resources/Actinos_6species.gbff"
    output:
        gff = "data/interim/prokka/{strains}/{strains}.gff",
        gbk = "data/interim/prokka/{strains}/{strains}.gbk",
        genus = temp("data/interim/prokka/{strains}/genus"),
        species = temp("data/interim/prokka/{strains}/species"),
    conda:
        "../envs/prokka.yaml"
    params:
        genus = "Streptomyces",
        increment = 10, 
        evalue = "1e-05"
    threads: 8
    shell:
        """
        head -1 {input.fna} | cut -d' ' -f2 > {output.genus}
        head -1 {input.fna} | cut -d' ' -f3 > {output.species}
        if [ `cat {output.genus}` == "Streptomyces" ]
        then
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --proteins {input.refseq} --prefix {wildcards.strains} --genus `cat {output.genus}` --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
        else
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --prefix {wildcards.strains} --genus `cat {output.genus}` --species `cat {output.species}` --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input}
        fi        
        """