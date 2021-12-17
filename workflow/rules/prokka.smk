rule copy_custom_fasta:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        cp {input} {output}
        """

rule extract_meta_prokka:
    input:
        fna = expand("data/interim/fasta/{strains}.fna", strains = STRAINS),
        samples_path = config["samples"],
        prokka_dir = "data/interim/prokka/",
        all_reports = expand("data/interim/assembly_report/{ncbi}.txt", ncbi = NCBI),
        assembly_report_path = "data/interim/assembly_report/",
    output:
        ncbi_meta_path = "data/processed/tables/df_ncbi_meta.csv",
        org_info = expand("data/interim/prokka/{strains}/organism_info.txt", strains = STRAINS),
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_organism_info.py {input.samples_path} {input.assembly_report_path} {input.prokka_dir} {output.ncbi_meta_path}
        """ 

if PROKKA_DB == []:
    rule prokka_default:
        input: 
            fna = "data/interim/fasta/{strains}.fna",
            org_info = "data/interim/prokka/{strains}/organism_info.txt"
        output:
            gff = "data/interim/prokka/{strains}/{strains}.gff",
            faa = "data/interim/prokka/{strains}/{strains}.faa",
            gbk = "data/interim/prokka/{strains}/{strains}.gbk",
        conda:
            "../envs/prokka.yaml"
        params:
            increment = 10, 
            evalue = "1e-05",
            rna_detection = "" # To use rnammer change value to --rnammer
        threads: 8
        log : "workflow/report/{strains}/prokka_run.log"
        shell:
            """
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --prefix {wildcards.strains} --genus "`cut -d "," -f 1 {input.org_info}`" --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} --evalue {params.evalue} {input.fna}
            cp data/interim/prokka/{wildcards.strains}/{wildcards.strains}.log {log}
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
            cat resources/prokka_db/gbk/*.gbff >> {output.refgbff}
            """

    rule prokka_custom:
        input: 
            fna = "data/interim/fasta/{strains}.fna",
            refgbff = "resources/prokka_db/reference.gbff",
            org_info = "data/interim/prokka/{strains}/organism_info.txt"
        output:
            gff = "data/interim/prokka/{strains}/{strains}.gff",
            faa = "data/interim/prokka/{strains}/{strains}.faa",
            gbk = "data/interim/prokka/{strains}/{strains}.gbk",
        conda:
            "../envs/prokka.yaml"
        params:
            increment = 10, 
            evalue = "1e-05",
            rna_detection = "" # To use rnammer change value to --rnammer
        threads: 8
        log : "workflow/report/{strains}/prokka_run.log"
        shell:
            """
            prokka --outdir data/interim/prokka/{wildcards.strains} --force --proteins {input.refgbff} --prefix {wildcards.strains} --genus "`cut -d "," -f 1 {input.org_info}`" --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} --evalue {params.evalue} {input.fna}
            cp data/interim/prokka/{wildcards.strains}/{wildcards.strains}.log {log}
            """

rule format_gbk:
    input: 
        gbk_prokka = "data/interim/prokka/{strains}/{strains}.gbk",
        gtdb = "data/processed/tables/df_gtdb_meta.csv",
    output:
        gbk_processed = "data/processed/genbank/{strains}.gbk",
    conda:
        "../envs/prokka.yaml"
    params:
        version = __version__
    shell:
        """
        python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} {params.version} {input.gtdb} {wildcards.strains} {output.gbk_processed}
        """
