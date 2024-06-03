if len(CUSTOM_GENBANK) > 0:
    rule copy_custom_genbank:
        input:
            gbk = lambda wildcards: DF_SAMPLES.loc[wildcards, "input_file"]
        output:
            gbk = "data/interim/processed-genbank/{custom_genbank}.gbk",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/prokka/copy_custom_fasta/copy_custom_fasta-{custom_genbank}.log"
        shell:
            """
            if [[ {input} == *.gb || {input} == *.gbk || {input} == *.genbank || {input} == *.gbff ]]
            then
                cp {input} {output} 2>> {log}
            else
                echo "ERROR: Wrong Extension:" {input} >> {log}
                exit 1
            fi
            """

    rule genbank_to_fna:
        input:
            gbk = "data/interim/processed-genbank/{custom_genbank}.gbk",
        output:
            fna = "data/interim/fasta/{custom_genbank}.fna",
        conda:
            "../envs/convert_genbank.yaml"
        log: "logs/convert_genbank/genbank_to_fna/genbank_to_fna-{custom_genbank}.log"
        shell:
            """
            any2fasta {input.gbk} > {output.fna} 2>> {log}
            """

    rule genbank_to_faa:
        input:
            gbk = "data/interim/processed-genbank/{custom_genbank}.gbk",
        output:
            faa = "data/interim/prokka/{custom_genbank}/{custom_genbank}.faa",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/convert_genbank/genbank_to_faa/genbank_to_faa-{custom_genbank}.log"
        shell:
            """
            python workflow/bgcflow/bgcflow/misc/create_aa.py {input.gbk} {output.faa} 2>> {log}
            """

    rule extract_meta_genbank:
        input:
            fna = "data/interim/fasta/{custom_genbank}.fna",
        output:
            org_info = "data/interim/prokka/{custom_genbank}/organism_info.txt",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/convert_genbank/extract_meta_prokka/extract_meta_prokka-{custom_genbank}.log"
        params:
            samples_path = bgcflow_util_dir / "samples.csv",
        shell:
            """
            python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.custom_genbank} \
                "{params.samples_path}" data/interim/assembly_report/ data/interim/prokka/ 2>> {log}
            """

    rule genbank_to_gff:
        input:
            gbk = "data/interim/processed-genbank/{custom_genbank}.gbk",
        output:
            gff = "data/interim/prokka/{custom_genbank}/{custom_genbank}.gff",
        conda:
            "../envs/convert_genbank.yaml"
        log: "logs/convert_genbank/genbank_to_gff/genbank_to_gff-{custom_genbank}.log"
        params:
            write_fasta = "True"
        shell:
            """
            python workflow/bgcflow/bgcflow/misc/convert_gbk_to_gff.py {input.gbk} {output.gff} {params.write_fasta} 2>> {log}
            """

    rule copy_converted_gbk:
        input:
            gbk = "data/interim/processed-genbank/{custom_genbank}.gbk",
        output:
            gbk = report("data/processed/{name}/genbank/{custom_genbank}.gbk", \
                caption="../report/file-genbank.rst", category="{name}", subcategory="Annotated Genbanks"),
            #tsv = "data/processed/{name}/genbank/{strains_fna}.tsv",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/custom_genbank/copy_converted_gbk/copy_converted_gbk_-{custom_genbank}-{name}.log"
        shell:
            """
            cp {input.gbk} {output.gbk} 2>> {log}
            """

    rule summarize_converted_gbk:
        input:
            gbk = "data/processed/{name}/genbank/{custom_genbank}.gbk"
        output:
            summary = "data/processed/{name}/genbank/{custom_genbank}.txt",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "logs/custom_genbank/summarize_converted_gbk/summarize_converted_gbk_-{custom_genbank}-{name}.log"
        shell:
            """
            python workflow/bgcflow/bgcflow/misc/summarize_gbk_txt.py {input.gbk} {output.summary} 2>> {log}
            """
