rule prepare_for_annotation:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    run:
        """
        cat {input} | sed "s/contig/{wildcards.custom}/" | sed "s/scaffold/{wildcards.custom}_scaf/" > {output}
        """

rule extract_meta_prokka:
    input:
        expand("data/interim/fasta/{strains}.fna", strains = STRAINS)
    output:
        expand("data/interim/prokka/{strains}/organism_info.txt", strains = STRAINS)
    run:
        # set up sample for default case with fasta files provided
        df_samples = pd.read_csv(config["samples"]).set_index("genome_id", drop=False)
        df_samples.index.names = ["genome_id"]

        for idx in df_samples.index:
            GENUS = df_samples.loc[idx, 'genus']
            SPECIES = df_samples.loc[idx, 'species']
            STRAIN_ID = df_samples.loc[idx, 'strain']
            org_info_path = os.path.join('data/interim/prokka/', idx, 'organism_info.txt')
            with open(org_info_path, 'w') as file_obj:
                file_obj.write(' '.join([GENUS,  SPECIES, STRAIN_ID]))

if PROKKA_DB == []:
    rule prokka_default:
        input: 
            fna = "data/interim/fasta/{strains}.fna",
            org_info = "data/interim/prokka/{strains}/organism_info.txt"
        output:
            gff = "data/interim/prokka/{strains}/{strains}.gff",
            gbk = "data/interim/prokka/{strains}/{strains}.gbk",
        conda:
            "../envs/prokka.yaml"
        params:
            increment = 10, 
            evalue = "1e-05"
        threads: 16
        shell:
            """
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --prefix {wildcards.strains} --genus `cut -d " " -f 1 {input.org_info}` --species `cut -d " " -f 2 {input.org_info}` --strain `cut -d " " -f 3 {input.org_info}` --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
            """
else:
    rule prokka_reference_download:
        output:
            gbff = "resources/prokka_db/gbk/{prokka_db}.gbff"
        conda:
            "../envs/prokka.yaml"
        shell:
            """
            ncbi-genome-download -s genbank -F genbank -A {wildcards.prokka_db} -o resources/prokka_db/download bacteria --verbose
            gunzip -c resources/prokka_db/download/genbank/bacteria/{wildcards.prokka_db}/*.gbff.gz > {output.gbff}
            rm -rf resources/prokka_db/download/genbank/bacteria/{wildcards.prokka_db}
            """

    rule prokka_db_setup:
        input:
            gbff = expand("resources/prokka_db/gbk/{prokka_db}.gbff", prokka_db = PROKKA_DB)
        output:
            refgbff = "resources/prokka_db/reference.gbff"
        conda:
            "../envs/prokka.yaml"
        shell:
            """
            cat resources/prokka_db/gbk/*.gbff > {output.refgbff}
            """

    rule prokka_custom:
        input: 
            fna = "data/interim/fasta/{strains}.fna",
            refgbff = "resources/prokka_db/reference.gbff",
            org_info = "data/interim/prokka/{strains}/organism_info.txt"
        output:
            gff = "data/interim/prokka/{strains}/{strains}.gff",
            gbk = "data/interim/prokka/{strains}/{strains}.gbk",
            gbk_processed = "data/processed/genbank/{strains}.gbk",
        conda:
            "../envs/prokka.yaml"
        params:
            increment = 10, 
            evalue = "1e-05"
        threads: 16
        shell:
            """
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --proteins {input.refgbff} --prefix {wildcards.strains} --genus `cut -d " " -f 1 {input.org_info}` --species `cut -d " " -f 2 {input.org_info}` --strain `cut -d " " -f 3 {input.org_info}` --cdsrnaolap --cpus {threads} --rnammer --increment {params.increment} --evalue {params.evalue} {input.fna}
            cp {output.gbk} {output.gbk_processed}
            """