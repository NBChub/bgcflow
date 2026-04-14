rule bigscape2:
    input:
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        antismash_dir="data/interim/bgcs/{name}/{version}/",
    output:
        index="data/interim/bigscape2/{name}_antismash_{version}/index.html",
    conda:
        "../envs/bigscape2.yaml"
    params:
        bigscape2_dir="data/interim/bigscape2/{name}_antismash_{version}/",
        label="{name}_antismash_{version}",
        pfam="resources/antismash_db/pfam/35.0/Pfam-A.hmm",
    log:
        "logs/bigscape2/{name}_antismash_{version}/bigscape2.log",
    threads: 64
    shell:
        """
        bigscape cluster -i {input.antismash_dir} -o {params.bigscape2_dir} -c {threads} --gcf-cutoffs 0.3,0.4,0.5 --include-singletons --label {params.label} --hybrids-off --mibig-version 4.0 -p {params.pfam} --verbose &>> {log}
        """


rule copy_bigscape2_zip:
    input:
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        index="data/interim/bigscape2/no_mibig_{name}_antismash_{version}/index.html",
    output:
        zip="data/processed/{name}/bigscape2/no_mibig_bigscape2_as{version}.zip",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape2/copy_bigscape2_zip/copy_bigscape2_zip-{name}-{version}.log",
    shell:
        """
        topdir=$PWD
        (cd data/interim/bigscape2 && zip -r $topdir/{output.zip} no_mibig_{wildcards.name}_antismash_{wildcards.version} \
            -x {wildcards.name}_antismash_{wildcards.version}/cache/**\* &>> $topdir/{log})
        """

rule bigscape2_to_cytoscape:
    input:
        index="data/interim/bigscape2/{name}_antismash_{version}/index.html",
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        mibig_bgc_table="resources/mibig/df_mibig_bgcs.csv/",
        as_dir="data/interim/bgcs/{name}/{version}",
        df_genomes_path="data/processed/{name}/tables/df_antismash_{version}_summary.csv",
    output:
        output_dir=directory(
            "data/processed/{name}/bigscape2/for_cytoscape_antismash_{version}"
        ),
        bgc_mapping="data/processed/{name}/bigscape2/{name}_bigscape2_as_{version}_mapping.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape2/bigscape2_to_cytoscape/bigscape2_to_cytoscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_bigscape2_to_cytoscape.py data/interim/bigscape2/{wildcards.name}_antismash_{wildcards.version}/ {input.as_dir} {input.df_genomes_path} {input.mibig_bgc_table} {output.output_dir} 2>> {log}
        python workflow/bgcflow/bgcflow/data/convert_region_to_bigscape2_region.py {input.bgc_mapping} {output.bgc_mapping} --overwrite 2>> {log}
        """


rule copy_bigscape2:
    input:
        index="data/interim/bigscape2/{name}_antismash_{version}/index.html",
        cytoscape="data/processed/{name}/bigscape2/for_cytoscape_antismash_{version}",
    output:
        main="data/processed/{name}/bigscape2/result_as{version}/index.html",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape2/copy_bigscape2_zip/copy_bigscape2-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/bigscape_copy.py {input.index} data/processed/{wildcards.name}/bigscape2/result_as{wildcards.version} 2>> {log}
        """
