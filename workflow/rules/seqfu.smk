rule seqfu_stats:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES)
    output:
        all_csv = report("data/processed/{name}/tables/df_seqfu_stats.csv", caption="../report/table-seqfu.rst", category="{name}", subcategory="Quality Control")
    conda:
        "../envs/seqfu.yaml"
    log: "workflow/report/logs/seqfu/seqfu-{name}.log"
    params:
        precision = 3,
    shell:
        """
        seqfu stats {input.fna} --csv -b --gc --precision 3 > {output.all_csv} 2>> {log}
        """