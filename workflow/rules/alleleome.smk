rule prepare_alleleome:
    input:
        roary="data/interim/roary/{name}"
    output:
        pangene_v2="data/processed/{name}/alleleome/pangene_v2.csv",
    log:
        "logs/prepare_alleleome/{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_reassign_pangene_categories.py \
            --data-dir {input.roary} \
            --output-file {output.pangene_v2} 2>> {log}
        """

rule prepare_alleleome_fasta:
    input:
        roary_path="data/interim/roary/{name}",
        gbk_folder="data/interim/processed-genbank",
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        fasta=directory("data/processed/{name}/pangenome_alignments")
    log:
        "logs/prepare_alleleome_fasta/{name}.log"
    conda:
        "../envs/alleleome.yaml"
    shell:
        """
        python workflow/scripts/alleleome_get_core_genes_fasta.py \
            --roary_path {input.roary_path} \
            --gbk_folder {input.gbk_folder} \
            --pangene_summary_path {input.pangene_summary_path} \
            --output_folder data/processed/{wildcards.name} 2>> {log}
        """

rule alleleome:
    input:
        pangenome_aligments_path="data/processed/{name}/pangenome_alignments",
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        alleleome_directory_path="data/processed/{name}/alleleome/pan_core_gene_syno_non_syno_df.csv"
    conda:
        "../envs/alleleome.yaml"
    log:
        "logs/alleleome/{name}.log"
    shell:
        """
        Alleleome --path1 {input.pangenome_aligments_path} \
            --path2 data/processed/{wildcards.name}/alleleome/ \
            --table {input.pangene_summary_path} \
            --log_to_terminal 2>> {log}
        """
