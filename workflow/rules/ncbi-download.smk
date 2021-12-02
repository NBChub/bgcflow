if NCBI == []:
    pass
else:
    rule ncbi_genome_download:
        output:
            gbff = "data/processed/genbank/{ncbi}.gbff",
            assembly_report = "data/interim/assembly_report/{ncbi}.txt"
        conda:
            "../envs/prokka.yaml"
        params:
            src = "genbank",
            fmt = "genbank"
        shell:
            """
            ncbi-genome-download -s {params.src} -F {params.fmt},assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose
            gunzip -c data/raw/ncbi/download/{params.src}/bacteria/{wildcards.ncbi}/*.gbff.gz > {output.gbff}
            cp data/raw/ncbi/download/{params.src}/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/{params.src}/bacteria/{wildcards.ncbi}
            """

    rule ncbi_metadata:
        input: 
            _all_ = expand("data/interim/assembly_report/{ncbi}.txt", ncbi = NCBI),
            assembly_report_path = "data/interim/assembly_report/",
        output:
            meta_out_path = "data/processed/tables/df_ncbi_meta.csv",
        conda:
            "../envs/prokka.yaml"
        script:
            "../bgcflow/bgcflow/data/make_ncbi_metadata.py"