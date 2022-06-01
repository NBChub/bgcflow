rule fastani:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
    output:
        fastani_infile = "data/interim/fastani/{name}/fastani_in.txt",
        fastani_out = "data/interim/fastani/{name}/fastani_out.tsv",
        fastani_matrix = "data/interim/fastani/{name}/fastani_out.tsv.matrix",
    conda:
        "../envs/fastani.yaml"
    threads: 32
    log: "workflow/report/logs/fastani/fastani-{name}.log"
    shell:
        """
        for fna in {input.fna} 
        do
            echo $fna >> {output.fastani_infile}
        done
        fastANI --ql {output.fastani_infile} --rl {output.fastani_infile} -t {threads} --matrix -o {output.fastani_out} 2>> {log} 
        """
