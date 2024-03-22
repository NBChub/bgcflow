rule prepare_pan_alleleome:
    input:
        roary="data/interim/roary/{name}"
    output:
        pangene_v2="data/processed/{name}/alleleome/pangene_v2.csv",
    log:
        "logs/prepare_pan_alleleome/{name}.log"
    conda:
        "../envs/pan_alleleome.yaml"
    shell:
        """
        python workflow/scripts/pan_alleleome_reassign_pangene_categories.py \
            --data-dir {input.roary} \
            --output-file {output.pangene_v2} 2>> {log}
        """

rule prepare_pan_alleleome_fasta:
    input:
        roary_path="data/interim/roary/{name}",
        gbk_folder=lambda wildcards: expand("data/interim/processed-genbank/{strains}.gbk",
            name=wildcards.name,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        fasta=directory("data/processed/{name}/pangenome_alignments")
    log:
        "logs/prepare_pan_alleleome_fasta/{name}.log"
    conda:
        "../envs/pan_alleleome.yaml"
    shell:
        """
        python workflow/scripts/pan_alleleome_get_pan_genes_fasta.py \
            --roary_path {input.roary_path} \
            --gbk_folder data/interim/processed-genbank \
            --pangene_summary_path {input.pangene_summary_path} \
            --output_folder data/processed/{wildcards.name} 2>> {log}
        """

rule pan_alleleome:
    input:
        pangenome_aligments_path="data/processed/{name}/pangenome_alignments",
        pangene_summary_path="data/processed/{name}/alleleome/pangene_v2.csv"
    output:
        alleleome_directory_path="data/processed/{name}/alleleome/pan_gene_syno_non_syno_df.csv"
    conda:
        "../envs/pan_alleleome.yaml"
    log:
        "logs/pan_alleleome/{name}.log"
    shell:
        """
        PAN_Alleleome --path1 {input.pangenome_aligments_path} \
            --path2 data/processed/{wildcards.name}/alleleome/ \
            --table {input.pangene_summary_path} \
            --log_to_terminal 2>> {log}
        """
