if NCBI == []:
    pass
else:
    rule ncbi_genome_download:
        output:
            fna = "data/interim/fasta/{ncbi}.fna",
            assembly_report = "data/interim/assembly_report/{ncbi}.txt"
        conda:
            "../envs/prokka.yaml"
        log: "workflow/report/logs/{ncbi}/ncbi_genome_download.log"
        shell:
            """
            ncbi-genome-download -s refseq -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose >> {log}
            gunzip -c data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*.fna.gz > {output.fna}
            cp data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/refseq/bacteria/{wildcards.ncbi}
            """