rule prepare_for_annotation:
    input:
        "results/genomes/{strains}/final_polish.fasta"
    output:
        "results/genomes/{strains}/prokkainput.fna" 
    shell:
        """
        cat {input} | sed "s/contig/{wildcards.strains}/" | sed "s/scaffold/{wildcards.strains}_scaf/" > {output}
        """

rule prokka:
    input: 
        "results/genomes/{strains}/prokkainput.fna",
        "results/genomes/{strains}/genus", 
        "results/genomes/{strains}/species",
        "resources/Actinos_6species.gbff"
    output:
        directory("results/genomes/{strains}/{strains}_prokka_actinoannotPFAM")
    conda:
        "../envs/prokka.yaml"
    threads: 12
    shell:
        """

        prokka --outdir {output} --proteins {input[3]} --prefix {wildcards.strains} --genus `cat {input[1]}` --species `cat {input[2]}` --strain {wildcards.strains} --cdsrnaolap --cpus {threads} --rnammer --increment 10 --evalue 1e-05 {input[0]}
        """
        
rule antismash:
    input: 
        "results/genomes/{strains}/{strains}_prokka_actinoannotPFAM",
        "resources/antismash_db.txt"
    output:
        directory("results/genomes/{strains}/{strains}_antiSMASH/"),
        "results/genomes/{strains}/{strains}_antiSMASH/{strains}.gbk"
    conda:
        "../envs/antismash.yaml"
    threads: 12
    shell:
        """
        antismash --genefinding-tool prodigal --output-dir {output[0]} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input[0]}/{wildcards.strains}.gbk -v
        """