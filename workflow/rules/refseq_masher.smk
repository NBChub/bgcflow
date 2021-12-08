rule refseq_masher:
    input: 
        fasta = "data/interim/fasta/{strains}.fna",
    output:
        masher_csv = "data/interim/refseq_masher/{strains}_masher.csv"
    conda:
        "../envs/refseq_masher.yaml"
    params:
        output_type = "csv",
        top_n_results = 10,
    shell:
        """
        refseq_masher matches --output-type {params.output_type} --top-n-results {params.top_n_results} --output {output.masher_csv} {input.fasta}
        """
