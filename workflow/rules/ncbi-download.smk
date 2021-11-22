if NCBI == []:
    pass
else:
    rule ncbi_genome_download:
        output:
            fna = "data/interim/fasta/{ncbi}.fna",
            assembly_report = "data/interim/assembly_report/{ncbi}.txt"
        conda:
            "../envs/prokka.yaml"
        shell:
            """
            ncbi-genome-download -s refseq -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose
            gunzip -c data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*.fna.gz > {output.fna}
            cp data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}
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
            "../src/data/make_ncbi_metadata.py"