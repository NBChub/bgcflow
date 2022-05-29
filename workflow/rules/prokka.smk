try:
    if Path(config["resources_path"]["RNAmmer"]).is_file():
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
except KeyError:
    prokka_params_rna = ""
    pass

rule copy_custom_fasta:
    input:
        "data/raw/fasta/{custom}.fna"
    output:
        "data/interim/fasta/{custom}.fna" 
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/copy_custom_fasta/copy_custom_fasta-{custom}.log"
    shell:
        """
        cp {input} {output} 2>> {log}
        """

rule prokka_db_setup:
    input:
        table = "resources/prokka_db/{prokka_db}.json"
    output:
        refgbff = "resources/prokka_db/{prokka_db}.gbff",
        ncbi_tempdir = temp(directory("resources/prokka_db/{prokka_db}")),
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/prokka_db_setup/prokka_db_setup-{prokka_db}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_prokka_db.py {wildcards.prokka_db} \
            {input.table} {output.refgbff} 2>> {log}
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
        python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.strains} \
            "{params.samples_path}" data/interim/assembly_report/ data/interim/prokka/ 2>> {log}
        """
try:
    if os.path.isfile(config["resources_path"]["pfam_for_prokka"]):
        prokka_use_pfam = f'--hmms {config["resources_path"]["pfam_for_prokka"]}'
    else:
        prokka_use_pfam = ""
except KeyError:
    prokka_use_pfam = ""

rule prokka:
    input: 
        fna = "data/interim/fasta/{strains}.fna",
        org_info = "data/interim/prokka/{strains}/organism_info.txt",
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, "file", DF_SAMPLES, PROKKA_DB_MAP)
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
        refgbff = lambda wildcards: get_prokka_refdb(wildcards, "params", DF_SAMPLES, PROKKA_DB_MAP),
        use_pfam = prokka_use_pfam
    threads: 4
    shell:
        """
        prokka --outdir data/interim/prokka/{wildcards.strains} --force \
            {params.refgbff} --prefix {wildcards.strains} --genus "`cut -d "," -f 1 {input.org_info}`" \
            --species "`cut -d "," -f 2 {input.org_info}`" --strain "`cut -d "," -f 3 {input.org_info}`" \
            --cdsrnaolap --cpus {threads} {params.rna_detection} --increment {params.increment} \
            {params.use_pfam} --evalue {params.evalue} {input.fna} &> {log}
        """

rule format_gbk:
    input: 
        gbk_prokka = "data/interim/prokka/{strains}/{strains}.gbk",
        gtdb_json = "data/interim/gtdb/{strains}.json",
    output:
        gbk_processed = "data/interim/processed-genbank/{strains}.gbk"
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        version = __version__,
    log: "workflow/report/logs/prokka/format_gbk/format_gbk-{strains}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/format_genbank_meta.py {input.gbk_prokka} \
            {params.version} {input.gtdb_json} {wildcards.strains} {output.gbk_processed} 2> {log}
        """

rule copy_prokka_gbk:
    input:
        gbk = "data/interim/processed-genbank/{strains}.gbk",
    output:
        gbk = report("data/processed/{name}/genbank/{strains}.gbk", \
            caption="../report/file-genbank.rst", category="{name}", subcategory="Annotated Genbanks")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/prokka/copy_prokka_gbk/copy_prokka_gbk_-{strains}-{name}.log"
    shell:
        """
        cp {input.gbk} {output.gbk} 2>> {log}
        """
