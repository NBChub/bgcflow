rule gtdb_prep:
    input:
        fasta = expand("data/interim/fasta/{strains}.fna", strains = STRAINS),
    output:
        taxonomy_placement = "data/interim/gtdb/placement_list.txt",
    log: "workflow/report/logs/gtdb_prep.log"
    conda: 
        "../envs/bgc_analytics.yaml"
    params:
        samples_path = SAMPLE_PATHS,
    shell: 
        """
        python workflow/bgcflow/bgcflow/data/gtdb_prep.py '{params.samples_path}' {output.taxonomy_placement} 2> {log}
        """

rule fetch_gtdb_taxonomy: #need to update by installing bigslice
    input: 
        taxonomy_placement = "data/interim/gtdb/placement_list.txt"
    output:
        taxonomy = temp("data/interim/gtdb/all_taxonomy_raw.tsv"),
    params:
        version = 'R202'    
    conda:
        "../envs/gtdb.yaml"
    log: "workflow/report/logs/fetch_gtdb_taxonomy.log"
    shell:
        """
        wget https://raw.githubusercontent.com/medema-group/bigslice/master/misc/assign_gtdb_taxonomy/fetch_taxonomy_from_api.py -P workflow/scripts/ -nc
        python workflow/scripts/fetch_taxonomy_from_api.py {input.taxonomy_placement} {output.taxonomy} --gtdb {params.version} 2> {log}
        """

rule fix_gtdb_taxonomy:
    input: 
        taxonomy_raw = "data/interim/gtdb/all_taxonomy_raw.tsv",
    output:
        taxonomy = "data/interim/gtdb/all_taxonomy.tsv",
        meta = "data/processed/tables/df_gtdb_meta.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/fix_gtdb_taxonomy.log"
    priority: 50
    params:
        samples_path = SAMPLE_PATHS,
    shell: 
        """
        python workflow/bgcflow/bgcflow/data/fix_gtdb_taxonomy.py '{params.samples_path}' {input.taxonomy_raw} {output.taxonomy} {output.meta} 2> {log}
        """