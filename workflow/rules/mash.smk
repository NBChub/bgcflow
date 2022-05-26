rule mash:
    input: 
        fna = lambda wildcards: get_fasta_inputs(wildcards.name),
    output:
        mash_infile = "data/interim/mash/{name}/mash_in.txt",
        triangle_dist = "data/interim/mash/{name}/triangle_distance_matrix.tsv"
    conda:
        "../envs/mash.yaml"
    threads: 32
    log: "workflow/report/logs/mash/mash-triangle-{name}.log"
    shell:
        """
        for fna in {input.fna} 
        do
            echo $fna >> {output.mash_infile}
        done
        (mash triangle -p {threads} -l {output.mash_infile} >> {output.triangle_dist}) 2>> {log}
        """
