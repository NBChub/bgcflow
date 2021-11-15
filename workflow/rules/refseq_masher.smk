rule refseq_masher:
    input: 
        fasta = expand("data/interim/fasta/{strains}.fna", strains=STRAINS)
    output:
        masher_csv = "data/processed/tables/refseq_masher.csv"
    conda:
        "../envs/refseq_masher.yaml"
    params:
        output_type = "csv",
        top_n_results = 10,
    shell:
        """
        refseq_masher matches --output-type {params.output_type} --top-n-results {params.top_n_results} --output {output.masher_csv} {input.fasta}
        """