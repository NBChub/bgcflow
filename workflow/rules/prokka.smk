if config["rules"]["rnammer"] == True:
    prokka_params_rna = "--rnammer"
    rule rnammer_setup:
        output: 
            "resources/rnammer_test.txt" 
        priority: 50
        conda:
            "../envs/prokka.yaml"
        log: "workflow/report/logs/prokka/rnammer_setup.log"
        shell:
            """
            ln -s $PWD/resources/RNAmmer/rnammer $CONDA_PREFIX/bin/rnammer 2>> {log}
            rnammer -S bac -m lsu,ssu,tsu -gff - resources/RNAmmer/example/ecoli.fsa >> {output}
            """
else:
    prokka_params_rna = ""
    pass

rule copy_custom_fasta:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/copy_custom_fasta-{custom}.log"
    shell:
        """
        cp {input} {output} 2>> {log}
        """

rule prokka_reference_download:
    output:
        gbff = temp("resources/prokka_db/gbk/{prokka_db}.gbff") # all projects ncbi accession 
    conda:
        "../envs/prokka.yaml"
    log: "workflow/report/logs/prokka/prokka_reference_download/prokka_reference_download-{prokka_db}.log"
    shell:
        """
        if [[ {wildcards.prokka_db} == GCF* ]]
        then
            source="refseq"
        elif [[ {wildcards.prokka_db} == GCA* ]]
        then
            source="genbank"
        else
            echo "accession must start with GCA or GCF" >> {log}
        fi
        ncbi-genome-download -s $source -F genbank -A {wildcards.prokka_db} -o resources/prokka_db/download bacteria --verbose >> {log}
        gunzip -c resources/prokka_db/download/$source/bacteria/{wildcards.prokka_db}/*.gbff.gz > {output.gbff}
        rm -rf resources/prokka_db/download/$source/bacteria/{wildcards.prokka_db}
        """

rule prokka_db_setup:
    input:
        gbff = get_prokka_db_accessions
    output:
        refgbff = "resources/prokka_db/reference_{name}.gbff"
    conda:
        "../envs/prokka.yaml"
    log: "workflow/report/logs/prokka/prokka_db_setup/prokka_db_setup-{name}.log"
    shell:
        """
        cat resources/prokka_db/gbk/*.gbff >> {output.refgbff}
        head {output.refgbff} >> {log}
        """

rule extract_meta_prokka:
    input:
        fna = "data/interim/fasta/{strains}.fna",
    output:
        org_info = "data/interim/prokka/{strains}/organism_info.txt",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/extract_meta_prokka/extract_meta_prokka-{strains}.log"
    params:
        samples_path = SAMPLE_PATHS,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.strains} "{params.samples_path}" data/interim/assembly_report/ data/interim/prokka/ 2>> {log}
        """

rule extract_ncbi_information:
    input:
        all_reports = expand("data/interim/prokka/{ncbi}/organism_info.txt", ncbi = NCBI),
        all_json = expand("data/interim/assembly_report/{ncbi}.json", ncbi = NCBI),
        assembly_report_path = "data/interim/assembly_report/",
    output:
        ncbi_meta_path = report("data/processed/tables/df_ncbi_meta.csv", caption="../report/table-ncbi_meta.rst", category="Genome Overview", subcategory="Metadata"),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/extract_ncbi_information.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/extract_ncbi_information.py {input.assembly_report_path} {output.ncbi_meta_path} 2>> {log}
        """

rule prokka:
    input: 
        fna = "data/interim/fasta/{strains}.fna",
        org_info = "data/interim/prokka/{strains}/organism_info.txt",
        refgbff = expand("resources/prokka_db/reference_{name}.gbff", name=PROJECT_IDS)
    output:
        gff = "data/interim/prokka/{strains}/{strains}.gff",
        faa = "data/interim/prokka/{strains}/{strains}.faa",
        gbk = "data/interim/prokka/{strains}/{strains}.gbk",
    conda: 
        "../envs/prokka.yaml"
    log: "workflow/report/logs/prokka/prokka/prokka-{strains}.log"
    params:
        increment = 10, 
        evalue = "1e-05",
        rna_detection = prokka_params_rna,
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, DF_SAMPLES)
    threads: 16
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.strains} --force {params.refgbff} --prefix {wildcards.strains} --genus "`cut -d "," -f 1 {input.org_info}`" --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} --evalue {params.evalue} {input.fna}
        cat data/interim/prokka/{wildcards.strains}/{wildcards.strains}.log > {log}
        """

rule format_gbk:
    input: 
        gbk_prokka = "data/interim/prokka/{strains}/{strains}.gbk",
        gtdb_json = "data/interim/gtdb/{strains}.json",
    output:
        gbk_processed = report("data/processed/genbank/{strains}.gbk", caption="../report/file-genbank.rst", category="Genome Overview", subcategory="Annotated Genbanks")
    conda:
        "../envs/prokka.yaml"
    params:
        version = __version__,
    log: "workflow/report/logs/prokka/format_gbk/format_gbk-{strains}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} {params.version} {input.gtdb_json} {wildcards.strains} {output.gbk_processed} 2> {log}
        """