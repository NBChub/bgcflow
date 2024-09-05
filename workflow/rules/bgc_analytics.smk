rule bgc_count:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.gbk",
    output:
        bgc_count="data/interim/antismash/{version}/{strains}_bgc_counts.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/bgc_counts/as_{version}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_bgc_counts.py {input.antismash} {wildcards.strains} {output.bgc_count} 2>> {log}
        """

rule antismash_overview:
    input:
        antismash="data/interim/antismash/{version}/{strains}/{strains}.json",
    output:
        bgc_table=temp("data/interim/antismash/{version}/{strains}_bgc_overview.json"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/bgc_overview/as_{version}_{strains}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_antismash_overview.py {input.antismash} {output.bgc_table} {wildcards.strains} 2>> {log}
        """

rule antismash_overview_gather:
    input:
        bgc_overview=lambda wildcards: expand("data/interim/antismash/{version}/{strains}_bgc_overview.json",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        mapping_dir = "data/interim/bgcs/{name}/{version}",
    output:
        df_bgc="data/processed/{name}/tables/df_regions_antismash_{version}.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/antismash_overview_gather-{version}-{name}.log",
    shell:
        """
        TMPDIR="data/interim/bgcs/{wildcards.name}/tmp/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_regions_antismash.txt"
        echo '{input.bgc_overview}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/get_antismash_overview_gather.py \
            $INPUT_JSON {input.mapping_dir} {output.df_bgc} 2>> {log}
        rm $INPUT_JSON
        """

rule copy_log_changes:
    input:
        mapping_dir="data/interim/bgcs/{name}/{version}",
    output:
        mapping_dir=directory("data/processed/{name}/log_changes/{version}")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/bgcs/downstream_bgc_prep/{name}/copy_log_changes-{version}.log",
    shell:
        """
        mkdir -p {output.mapping_dir} 2>> {log}
        cp {input.mapping_dir}/*/*-change_log.json {output.mapping_dir}/. 2>> {log}
        """

rule antismash_summary:
    input:
        mapping_dir="data/processed/{name}/log_changes/{version}",
        bgc_count=lambda wildcards: expand(
            "data/interim/antismash/{version}/{strains}_bgc_counts.json",
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        as_dir=lambda wildcards: expand(
            "data/processed/{name}/antismash/{version}/{strains}",
            name=wildcards.name,
            version=wildcards.version,
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
        bgc_overview="data/processed/{name}/tables/df_regions_antismash_{version}.csv",
    output:
        df_antismash_summary=report(
            "data/processed/{name}/tables/df_antismash_{version}_summary.csv",
            caption="../report/table-antismash.rst",
            category="{name}",
            subcategory="AntiSMASH Summary Table",
        )
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/antismash_summary-{version}-{name}.log",
    params:
        df_samples=lambda wildcards: PEP_PROJECTS[wildcards.name].config["sample_table"],
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}/{wildcards.version}"
        mkdir -p $TMPDIR
        INPUT_JSON="$TMPDIR/df_bgc_counts.txt"
        echo '{input.bgc_count}' > $INPUT_JSON
        python workflow/bgcflow/bgcflow/data/make_genome_dataset.py $INPUT_JSON '{params.df_samples}' {output.df_antismash_summary} 2>> {log}
        rm $INPUT_JSON
        """

rule write_dependency_versions:
    input:
        "data/processed/{name}/metadata/project_metadata.json",
    output:
        "data/processed/{name}/metadata/dependency_versions.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/{name}/write_dependency_versions.log",
    params:
        antismash_version=antismash_major_version
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_dependencies.py {output} {params.antismash_version} 2> {log}
        """

rule get_project_metadata:
    output:
        "data/processed/{name}/metadata/project_metadata.json",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bgc_analytics/get_{name}_metadata.log",
    params:
        bgcflow_version=__version__,
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_project_metadata.py {wildcards.name} {output} {params.bgcflow_version} 2> {log}
        """

rule get_mibig_table:
    output:
        mibig_json_folder=directory("resources/mibig/json/"),
        mibig_bgc_table="resources/mibig/df_mibig_bgcs.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/get_mibig_table.log",
    params:
        mibig_version="3.1",
    shell:
        """
        # Define the MIBiG tar.gz file path
        TAR_GZ_FILE="resources/mibig_json_{params.mibig_version}.tar.gz"

        # Check if the tar.gz file exists and remove it if it does
        if [ -f "$TAR_GZ_FILE" ]; then
            rm "$TAR_GZ_FILE"
        fi

        # Download the tar.gz file
        (cd resources && wget https://dl.secondarymetabolites.org/mibig/mibig_json_{params.mibig_version}.tar.gz -nc) &>> {log}

        # Extract the tar.gz file and organize the contents
        (cd resources && tar -xvf mibig_json_{params.mibig_version}.tar.gz && mkdir -p mibig && mv mibig_json_{params.mibig_version}/ mibig/json && rm mibig_json_{params.mibig_version}.tar.gz) &>> {log}

        # Run the Python script to convert into table
        python workflow/bgcflow/bgcflow/data/get_mibig_data.py {output.mibig_json_folder} {output.mibig_bgc_table} 2>> {log}
        """
