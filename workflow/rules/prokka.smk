rule prepare_for_annotation:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    shell:
        """
        cat {input} | sed "s/contig/{wildcards.custom}/" | sed "s/scaffold/{wildcards.custom}_scaf/" > {output}
        """

rule prokka_genusdb_setup:
    input:
        lambda wildcards: PROKKA_DB_FILE[wildcards.prokka_genus]
    output:
        folder = directory("resources/prokka-db/{prokka_genus}"),
        report = "resources/prokka-db/{prokka_genus}_db.txt"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s genbank -F genbank -A {input} -o {output.folder} bacteria --verbose
        while read ln; do gunzip -c {output.folder}/genbank/bacteria/$ln/*.gbff.gz > {output.folder}/genbank/$ln.gbff; done < {input}
        prokka-genbank_to_fasta_db {output.folder}/genbank/*.gbff > {output.folder}/{wildcards.prokka_genus}.faa
        cd-hit -i {output.folder}/{wildcards.prokka_genus}.faa -o {output.folder}/{wildcards.prokka_genus} -T 0 -M 0 -g 1 -s 0.8 -c 0.9
        rm -fv {output.folder}/{wildcards.prokka_genus}.faa {output.folder}/{wildcards.prokka_genus}.clstr
        makeblastdb -dbtype prot -in {output.folder}/{wildcards.prokka_genus}
        mv {output.folder}/{wildcards.prokka_genus}.p* $CONDA_PREFIX/db/genus/
        echo $CONDA_PREFIX/db/genus/ > resources/prokka-db/{wildcards.prokka_genus}_location.txt
        ls $CONDA_PREFIX/db/genus/{wildcards.prokka_genus} >> resources/prokka-db/{wildcards.prokka_genus}_db.txt
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
<<<<<<< HEAD
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --usegenus --prefix {wildcards.strains} --genus `cat {output.genus}` --species `cat {output.species}` --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input}
=======
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --prefix {wildcards.strains} --genus `cat {output.genus}` --species `cat {output.species}` --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input}
>>>>>>> f785353ad46fc8cbca06ea04e5da837b477243ed
        fi        
        """