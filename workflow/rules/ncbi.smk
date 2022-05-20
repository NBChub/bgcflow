if NCBI == []:
    pass
else:
    rule ncbi_genome_download:
        output:
            fna = "data/interim/fasta/{ncbi}.fna",
            assembly_report = "data/interim/assembly_report/{ncbi}.txt",
            json_report = "data/interim/assembly_report/{ncbi}.json",
        conda:
            "../envs/bgc_analytics.yaml"
        log: "workflow/report/logs/ncbi/ncbi_genome_download/ncbi_genome_download_{ncbi}.log"
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

    rule extract_ncbi_information:
        input:
            all_json = lambda wildcards: get_ncbi_assembly_inputs(wildcards.name)
        output:
            ncbi_meta_path = report("data/processed/{name}/tables/df_ncbi_meta.csv", \
                caption="../report/table-ncbi_meta.rst", \
                category="Genome Overview", subcategory="Metadata"),
        conda:
            "../envs/bgc_analytics.yaml"
        log: "workflow/report/logs/ncbi/extract_ncbi_information/extract_ncbi_information-{name}.log"
        shell:
            """
            echo 1 2>> {log}
            python workflow/bgcflow/bgcflow/data/extract_ncbi_information.py \
                '{input.all_json}' {output.ncbi_meta_path} 2>> {log}
        """