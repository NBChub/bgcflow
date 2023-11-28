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
            if [[ {input} == *.gb || {input} == *.gbk || {input} == *.genbank ]]
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
        log: "logs/genbank_to_fna/genbank_to_fna/genbank_to_fna-{custom_genbank}.log"
        shell:
            """
            any2fasta {input.gbk} > {output.fna} 2>> {log}
            """
