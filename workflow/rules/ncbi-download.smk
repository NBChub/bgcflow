rule ncbi_genome_download:
    output:
        fna = "data/interim/fasta/{ncbi}.fna",
        assembly_report = "data/interim/assembly_report/{ncbi}.txt"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s refseq -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose
        gunzip -c data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*.fna.gz > {output.fna}
        cp data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
        rm -rf data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}
        """

rule ncbi_metadata:
    input: 
        _all_ = expand("data/interim/assembly_report/{ncbi}.txt", ncbi = NCBI),
        assembly_report_path = "data/interim/assembly_report/",
    output:
        meta_out_path = "data/processed/tables/df_ncbi_meta.csv",
        genus = expand("data/interim/assembly_report/{ncbi}.genus", ncbi = NCBI),
        species = expand("data/interim/assembly_report/{ncbi}.species", ncbi = NCBI),
        strain = expand("data/interim/assembly_report/{ncbi}.strain", ncbi = NCBI),
    conda:
        "../envs/prokka.yaml"
    script:
        "../src/data/make_ncbi_metadata.py"

rule prokka_ncbi:
    input: 
        fna = "data/interim/fasta/{ncbi}.fna",
        meta_out_path = "data/processed/tables/df_ncbi_meta.csv",
        refseq = "resources/Actinos_6species.gbff",
        genus = "data/interim/assembly_report/{ncbi}.genus",
        species = "data/interim/assembly_report/{ncbi}.species",
        strain = "data/interim/assembly_report/{ncbi}.strain",
    output:
        gff = "data/interim/prokka/{ncbi}/{ncbi}.gff",
        gbk = "data/interim/prokka/{ncbi}/{ncbi}.gbk",
    conda:
        "../envs/prokka.yaml"
    params:
        increment = 10, 
        evalue = "1e-05"
    threads: 12
    shell:
        """
        if [ `cat {input.genus}` == "Streptomyces" ]
        then
            prokka --outdir data/interim/prokka/{wildcards.ncbi} --force --prefix {wildcards.ncbi} --genus `cat {input.genus}` --species `cat {input.species}` --strain `cat {input.strain}` --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} --proteins {input.refseq} {input.fna}
        else
            prokka --outdir data/interim/prokka/{wildcards.ncbi} --force --prefix {wildcards.ncbi} --genus `cat {input.genus}` --species `cat {input.species}` --strain `cat {input.strain}` --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
        fi
        rm {input.genus} {input.species} {input.strain}
        """