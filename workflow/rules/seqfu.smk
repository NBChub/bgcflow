rule seqfu_stats:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES)
    output:
        all_csv = report("data/processed/{name}/tables/df_seqfu_stats.csv", caption="../report/table-seqfu.rst", category="Quality Control")
    conda:
        "../envs/seqfu.yaml"
    log: "workflow/report/logs/seqfu/seqfu-{name}.log"
    shell:
        """
        seqfu stats {input.fna} --csv -b > {output.all_csv} 2>> {log}
        """