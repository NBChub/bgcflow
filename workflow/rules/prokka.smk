if config["rules"]["rnammer"] == True:
    prokka_params_rna = "--rnammer"
    rule rnammer_setup:
        output: 
            "resources/rnammer_test.txt" 
        priority: 50
        conda:
            "../envs/prokka.yaml"
        log: "workflow/report/logs/rnammer_setup.log"
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
    log: "workflow/report/logs/{custom}/prokka_custom_copy.log"
    shell:
        """
        cp {input} {output} 2>> {log}
        """

rule prokka_reference_download:
    output:
        gbff = temp("resources/prokka_db/gbk/{prokka_db}.gbff") # all projects ncbi accession 
    conda:
        "../envs/prokka.yaml"
    log: "workflow/report/logs/prokka_db/{prokka_db}.log"
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
    log: "workflow/report/logs/prokka_db/prokka_db_{name}.log"
    shell:
        """
        cat resources/prokka_db/gbk/*.gbff >> {output.refgbff}
        head {output.refgbff} >> {log}
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
    log: "workflow/report/logs/{strains}/prokka_run.log"
    params:
        increment = 10, 
        evalue = "1e-05",
        rna_detection = prokka_params_rna,
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, DF_SAMPLES)
    threads: 8
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.strains} --force {params.refgbff} --prefix {wildcards.strains} --genus "`cut -d "," -f 1 {input.org_info}`" --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} --evalue {params.evalue} {input.fna}
        cat data/interim/prokka/{wildcards.strains}/{wildcards.strains}.log > {log}
        """

rule extract_meta_prokka:
    input:
        fna = expand("data/interim/fasta/{strains}.fna", strains = STRAINS),
        prokka_dir = "data/interim/prokka/",
        all_reports = expand("data/interim/assembly_report/{ncbi}.txt", ncbi = NCBI),
        assembly_report_path = "data/interim/assembly_report/",
    output:
        ncbi_meta_path = "data/processed/tables/df_ncbi_meta.csv",
        org_info = expand("data/interim/prokka/{strains}/organism_info.txt", strains = STRAINS),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka_meta.log"
    params:
        samples_path = SAMPLE_PATHS,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_organism_info.py "{params.samples_path}" {input.assembly_report_path} {input.prokka_dir} {output.ncbi_meta_path} 2>> {log}
        """ 

rule format_gbk:
    input: 
        gbk_prokka = "data/interim/prokka/{strains}/{strains}.gbk",
    output:
        gbk_processed = "data/processed/genbank/{strains}.gbk",
    conda:
        "../envs/prokka.yaml"
    params:
        version = __version__,
        gtdb = "data/processed/tables/df_gtdb_meta.csv",
    log: "workflow/report/logs/{strains}/prokka_format_gbk.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} {params.version} {params.gtdb} {wildcards.strains} {output.gbk_processed} 2> {log}
        """