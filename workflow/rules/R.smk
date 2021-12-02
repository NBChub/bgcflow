rule mag_stats:
    input: 
        gtdb_ar = "workflow/Rscripts/gtdbtk.ar122.summary.tsv",
        gtdb_bac = "workflow/Rscripts/gtdbtk.bac120.summary.tsv"
    output:
        "data/workflow/Rscripts/mag_stats.csv"
    conda:
        "../envs/R.yaml"
    shell:
        """
        Rscript workflow/Rscripts/mag_stats.r -o data/workflow/Rscripts/ -a {input.gtdb_ar} -b{gtdb_bac}
        """