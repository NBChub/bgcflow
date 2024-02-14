if "taxonomic_mode" not in config.keys():
    config["taxonomic_mode"] = "bacteria"

if NCBI == []:
    pass
else:
    if len(NCBI_FASTA) > 0:
        rule ncbi_genome_download_fasta:
            output:
                fna="data/interim/fasta/{ncbi_fasta}.fna",
                assembly_report="data/interim/assembly_report/{ncbi_fasta}.txt",
                json_report="data/interim/assembly_report/{ncbi_fasta}.json",
            conda:
                "../envs/bgc_analytics.yaml"
            log:
                "logs/ncbi/ncbi_genome_download/ncbi_genome_download_{ncbi_fasta}.log",
            params:
                ncbi_groups=config["taxonomic_mode"],
                file_format="fasta",
                extension="fna",
                retries=3,
            shell:
                """
                set -e
                if [[ {wildcards.ncbi_fasta} == GCF* ]]
                then
                    source="refseq"
                elif [[ {wildcards.ncbi_fasta} == GCA* ]]
                then
                    source="genbank"
                else
                    echo "accession must start with GCA or GCF" >> {log}
                fi
                ncbi-genome-download -s $source -F {params.file_format},assembly-report -A {wildcards.ncbi_fasta} -o data/raw/ncbi/download -P -N --verbose -d -r {params.retries} {params.ncbi_groups} 2>> {log}
                gunzip -c data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_fasta}/*.{params.extension}.gz > {output.fna}
                cp data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_fasta}/*report.txt {output.assembly_report}
                rm -rf data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_fasta}
                python workflow/bgcflow/bgcflow/data/get_assembly_information.py {output.assembly_report} {output.json_report} {wildcards.ncbi_fasta} 2>> {log}
                """

    elif len(NCBI_GENBANK) > 0:
        rule ncbi_genome_download_gbk:
            output:
                fna="data/interim/fasta/{ncbi_genbank}.fna",
                gff="data/interim/prokka/{ncbi_genbank}/{ncbi_genbank}.gff",
                faa="data/interim/prokka/{ncbi_genbank}/{ncbi_genbank}.faa",
                gbk="data/interim/processed-genbank/{ncbi_genbank}.gbk",
                assembly_report="data/interim/assembly_report/{ncbi_genbank}.txt",
                json_report="data/interim/assembly_report/{ncbi_genbank}.json",
            conda:
                "../envs/bgc_analytics.yaml"
            log:
                "logs/ncbi/ncbi_genome_download/ncbi_genome_download_{ncbi_genbank}.log",
            params:
                ncbi_groups=config["taxonomic_mode"],
                file_format="genbank,fasta,gff,translated-cds",
                retries=3,
            shell:
                """
                set -e
                if [[ {wildcards.ncbi_genbank} == GCF* ]]
                then
                    source="refseq"
                elif [[ {wildcards.ncbi_genbank} == GCA* ]]
                then
                    source="genbank"
                else
                    echo "accession must start with GCA or GCF" >> {log}
                fi
                ncbi-genome-download -s $source -F {params.file_format},assembly-report -A {wildcards.ncbi_genbank} -o data/raw/ncbi/download -P -N --verbose -d -r {params.retries} {params.ncbi_groups} 2>> {log}
                gunzip -c data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}/*.gbff.gz > {output.gbk}
                gunzip -c data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}/*.fna.gz > {output.fna}
                gunzip -c data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}/*.gff.gz > {output.gff}
                gunzip -c data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}/*.faa.gz > {output.faa}
                cp data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}/*report.txt {output.assembly_report}
                rm -rf data/raw/ncbi/download/$source/{params.ncbi_groups}/{wildcards.ncbi_genbank}
                python workflow/bgcflow/bgcflow/data/get_assembly_information.py {output.assembly_report} {output.json_report} {wildcards.ncbi_genbank} 2>> {log}
                """

    rule extract_ncbi_information:
        input:
            all_json=lambda wildcards: get_ncbi_assembly_inputs(
                wildcards.name, DF_SAMPLES
            ),
        output:
            ncbi_meta_path=report(
                "data/processed/{name}/tables/df_ncbi_meta.csv",
                caption="../report/table-ncbi_meta.rst",
                category="{name}",
                subcategory="NCBI Genome Overview",
            ),
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/ncbi/extract_ncbi_information/extract_ncbi_information-{name}.log",
        shell:
            """
            TMPDIR="data/interim/tmp/{wildcards.name}"
            mkdir -p $TMPDIR
            INPUT_JSON="$TMPDIR/df_ncbi_meta_input.txt"
            echo '{input.all_json}' > $INPUT_JSON
            python workflow/bgcflow/bgcflow/data/extract_ncbi_information.py \
                $INPUT_JSON {output.ncbi_meta_path} 2>> {log}
            rm $INPUT_JSON
            """

    rule download_patric_tables:
        output:
            patric_genome_summary="resources/patric_meta/genome_summary",
            patric_genome_metadata="resources/patric_meta/genome_metadata",
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/patric/download_patric_tables.log",
        shell:
            """
            wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_summary -O {output.patric_genome_summary} 2>> {log}
            wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata -O {output.patric_genome_metadata} 2>> {log}
            """

    rule extract_patric_meta:
        input:
            ncbi_meta_path="data/processed/{name}/tables/df_ncbi_meta.csv",
            patric_genome_summary="resources/patric_meta/genome_summary",
            patric_genome_metadata="resources/patric_meta/genome_metadata",
        output:
            patric_meta_path=report(
                "data/processed/{name}/tables/df_patric_meta.csv",
                caption="../report/table-patric_meta.rst",
                category="{name}",
                subcategory="Patric Genome Overview",
            ),
        conda:
            "../envs/bgc_analytics.yaml"
        log:
            "logs/patric/extract_patric_information/extract_patric_information-{name}.log",
        shell:
            """
            python workflow/bgcflow/bgcflow/data/extract_patric_meta.py \
                {input.ncbi_meta_path} {input.patric_genome_summary} {input.patric_genome_metadata} {output.patric_meta_path} 2>> {log}
            """
