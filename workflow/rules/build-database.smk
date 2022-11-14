rule antismash_json_extract:
    input:
        json = "data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        cdss = "data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
        dna_sequences = "data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
        regions = "data/interim/database/as_{version}/{strains}/{strains}_regions.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/scatter/as_{version}_json_extract_{strains}.log"
    params:
        outdir = "data/interim/database/as_{version}/{strains}",
    shell:  
        """
        python workflow/bgcflow/bgcflow/database/bgc_meta.py {input.json} {params.outdir} 2>> {log}
        """

rule build_dna_sequences_table:
    input:
        dna_sequences = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_dna_sequences.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        dna_sequences = "data/processed/{name}/tables/df_dna_sequences_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_dna_sequences_gather_{name}.log"
    shell:  
        """
        python workflow/bgcflow/bgcflow/database/create_cdss_table.py '{input.dna_sequences}' {output.dna_sequences} 2>> {log}
        """

rule build_regions_table:
    input:
        regions = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_regions.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        regions = "data/processed/{name}/tables/df_regions_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_regions_gather_{name}.log"
    shell:  
        """
        python workflow/bgcflow/bgcflow/database/create_cdss_table.py '{input.regions}' {output.regions} 2>> {log}
        """

rule build_cdss_table:
    input:
        cdss = lambda wildcards: expand("data/interim/database/as_{version}/{strains}/{strains}_cdss.json",
                                         version=wildcards.version,
                                         strains=[s for s in list(PEP_PROJECTS[wildcards.name].sample_table.index)]),
    output:
        cdss = "data/processed/{name}/tables/df_cdss_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/gather/as_{version}_cdss_gather_{name}.log"
    shell:  
        """
        python workflow/bgcflow/bgcflow/database/create_cdss_table.py '{input.cdss}' {output.cdss} 2>> {log}
        """

### draft rule to build duckdb database
rule build_database:
    input:
        cdss = "data/processed/{name}/tables/df_cdss_{version}.csv",
        regions = "data/processed/{name}/tables/df_regions_{version}.csv",
        dna_sequences = "data/processed/{name}/tables/df_dna_sequences_{version}.csv",
    output:
        log = "data/processed/{name}/tables/database_{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/database/report/database_{version}_{name}.log"
    shell:  
        """
        echo {input.cdss} >> {output.log}
        echo {input.regions} >> {output.log}
        echo {input.dna_sequences} >> {output.log}
        """