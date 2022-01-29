if NCBI == []:
    pass
else:
    rule ncbi_genome_download:
        output:
            fna = "data/interim/fasta/{ncbi}.fna",
            assembly_report = "data/interim/assembly_report/{ncbi}.txt",
            json_report = "data/interim/assembly_report/{ncbi}.json",
        conda:
            "../envs/ncbi_utilities.yaml"
        log: "workflow/report/logs/ncbi-download/ncbi_genome_download_{ncbi}.log"
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
            ncbi-genome-download -s $source -F fasta,assembly-report -A {wildcards.ncbi} -o data/raw/ncbi/download -P -N --verbose bacteria 2>> {log}
            gunzip -c data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}/*.fna.gz > {output.fna}
            cp data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}/*report.txt {output.assembly_report}
            rm -rf data/raw/ncbi/download/$source/bacteria/{wildcards.ncbi}
            python workflow/bgcflow/bgcflow/data/get_assembly_information.py {output.assembly_report} {output.json_report} {wildcards.ncbi} 2>> {log}
            """