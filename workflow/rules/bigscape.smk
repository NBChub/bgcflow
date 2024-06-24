rule install_bigscape:
    output:
        bigscape=directory("resources/BiG-SCAPE/"),
    conda:
        "../envs/bigscape.yaml"
    log:
        "logs/bigscape/bigscape-install_bigscape.log",
    params:
        release="1.1.9",
    shell:
        """
        (cd resources && wget https://github.com/medema-group/BiG-SCAPE/archive/refs/tags/v{params.release}.zip) &>> {log}
        (cd resources && unzip -o v{params.release}.zip && mv BiG-SCAPE-{params.release}/ BiG-SCAPE/ && rm v{params.release}.zip) &>> {log}
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz) &>> {log}
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm) &>> {log}
        """


rule bigscape:
    input:
        bigscape="resources/BiG-SCAPE",
        bgc_mapping=rules.downstream_bgc_prep.output.bgc_mapping,
        antismash_dir=rules.downstream_bgc_prep.output.outdir,
    output:
        index="data/interim/bigscape/{name}_antismash_{version}/index.html",
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir="data/interim/bigscape/{name}_antismash_{version}/",
        label="{name}_antismash_{version}",
    log:
        "logs/bigscape/{name}_antismash_{version}/bigscape.log",
    threads: 32
    shell:
        """
        python {input.bigscape}/bigscape.py -i {input.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose &>> {log}
        """


## Unused rule, might be deleted
rule bigscape_no_mibig:
    input:
        bigscape="resources/BiG-SCAPE",
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        antismash_dir="data/interim/bgcs/{name}/{version}/",
    output:
        index="data/interim/bigscape/no_mibig_{name}_antismash_{version}/index.html",
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir="data/interim/bigscape/no_mibig_{name}_antismash_{version}/",
        label="{name}_antismash_{version}",
    log:
        "logs/bigscape_no_mibig/{name}_antismash_{version}/bigscape.log",
    threads: 32
    shell:
        """
        python {input.bigscape}/bigscape.py -i {input.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --verbose --no_classify --mix &>> {log}
        """

rule copy_bigscape_zip:
    input:
        bgc_mapping=rules.downstream_bgc_prep.output.bgc_mapping,
        index=rules.bigscape.output.index,
    output:
        zip="data/processed/{name}/bigscape/no_mibig_bigscape_as{version}.zip",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/copy_bigscape_zip/copy_bigscape_zip-{name}-{version}.log",
    shell:
        """
        topdir=$PWD
        (cd data/interim/bigscape && zip -r $topdir/{output.zip} no_mibig_{wildcards.name}_antismash_{wildcards.version} \
            -x {wildcards.name}_antismash_{wildcards.version}/cache/**\* &>> $topdir/{log})
        """


rule bigscape_to_cytoscape:
    input:
        index=rules.bigscape.output.index,
        bgc_mapping=rules.downstream_bgc_prep.output.bgc_mapping,
        mibig_bgc_table="resources/mibig/df_mibig_bgcs.csv",
        as_dir="data/interim/bgcs/{name}/{version}",
        df_genomes_path=rules.antismash_summary.output.df_antismash_summary,
    output:
        output_dir=temp(directory(
            "data/interim/bigscape/for_cytoscape/{name}/{version}"
        )),
        bgc_mapping="data/processed/{name}/bigscape/{name}_bigscape_as_{version}_mapping.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/bigscape_to_cytoscape/bigscape_to_cytoscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_bigscape_to_cytoscape.py data/interim/bigscape/{wildcards.name}_antismash_{wildcards.version}/ {input.as_dir} {input.df_genomes_path} {input.mibig_bgc_table} {output.output_dir} 2>> {log}
        cp {input.bgc_mapping} {output.bgc_mapping}
        """

rule check_missing_bigscape:
    input:
        processed_bigscape=rules.bigscape_to_cytoscape.output.output_dir,
        antismash_region_table=rules.antismash_overview_gather.output.df_bgc,
        bigscape_log=rules.bigscape.log,
    output:
        cytoscape=directory("data/processed/{name}/bigscape/for_cytoscape_antismash_{version}"),
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/check_missing_bigscape/check_missing_bigscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/check_missing_bigscape.py {input.processed_bigscape} \
            {input.antismash_region_table} \
            {input.bigscape_log} \
            {output.cytoscape} 2>> {log}
        """

rule copy_bigscape:
    input:
        index=rules.bigscape.output.index,
        cytoscape=rules.check_missing_bigscape.output.cytoscape,
    output:
        main="data/processed/{name}/bigscape/result_as{version}/index.html",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/bigscape/copy_bigscape_zip/copy_bigscape-{name}-{version}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/bigscape_copy.py {input.index} data/processed/{wildcards.name}/bigscape/result_as{wildcards.version} 2>> {log}
        """
