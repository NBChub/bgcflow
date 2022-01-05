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
            if [[ {wildcards.ncbi} == GCF* ]]
            then
                source="refseq"
            elif [[ {wildcards.ncbi} == GCA* ]]
            then
                source="genbank"
            else
                echo "accession must start with GCA or GCF" >> {log}
            fi
            ncbi-genome-download -s $source -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download bacteria --verbose >> {log}
            gunzip -c data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}/*.fna.gz > {output.fna}
            cp data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}
            """