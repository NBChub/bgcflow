rule fastani:
    input:
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fastani_infile="data/interim/fastani/{name}/fastani_in.txt",
        fastani_out="data/interim/fastani/{name}/fastani_out.tsv",
        fastani_matrix="data/interim/fastani/{name}/fastani_out.tsv.matrix",
    conda:
        "../envs/fastani.yaml"
    threads: 32
    log:
        "workflow/report/logs/fastani/fastani-{name}.log",
    shell:
        """
        for fna in {input.fna}
        do
            echo $fna >> {output.fastani_infile}
        done
        fastANI --ql {output.fastani_infile} --rl {output.fastani_infile} -t {threads} --matrix -o {output.fastani_out} 2>> {log}
        """


rule fastani_convert:
    input:
        fastani_matrix="data/interim/fastani/{name}/fastani_out.tsv.matrix",
    output:
        df_fastani="data/processed/{name}/fastani/df_fastani.csv",
    log:
        "workflow/report/logs/fastani/fastani-convert-{name}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/convert_triangular_matrix.py {input.fastani_matrix} {output.df_fastani} 2>> {log}
        """
