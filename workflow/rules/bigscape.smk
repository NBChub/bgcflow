rule install_bigscape:
    output:
        bigscape = directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/bigscape/bigscape-install_bigscape.log"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip) 2>> {log}
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/ && rm BiG-SCAPE-master.zip) 2>> {log}
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz) 2>> {log}
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm) 2>> {log}
        """

rule bigscape:
    input: 
        bigscape = "resources/BiG-SCAPE",
        bgc_mapping = "data/interim/bgcs/{name}/{name}_antismash_{version}.csv"
    output:
        index = "data/interim/bigscape/{name}_antismash_{version}/index.html"
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir = "data/interim/bigscape/{name}_antismash_{version}/",
        label = "{name}_antismash_{version}",
        antismash_dir = "data/interim/bgcs/{name}/{version}/"
    log: "workflow/report/logs/bigscape/{name}_antismash_{version}/bigscape.log"
    threads: 32
    shell:
        """
        python {input.bigscape}/bigscape.py -i {params.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}
        """

rule copy_bigscape_zip:
    input: 
        bgc_mapping = "data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
        index = "data/interim/bigscape/{name}_antismash_{version}/index.html"
    output:
        zip = "data/processed/{name}/bigscape/{name}_bigscape_as{version}.zip"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/bigscape/copy_bigscape_zip/copy_bigscape_zip-{name}-{version}.log"
    shell:
        """
        topdir=$PWD
        (cd data/interim/bigscape && zip -r $topdir/{output.zip} {wildcards.name}_antismash_{wildcards.version}.csv \
            {wildcards.name}_antismash_{wildcards.version} -x {wildcards.name}_antismash_{wildcards.version}/cache/**\* &>> $topdir/{log})
        """

rule bigscape_to_cytoscape:
    input:
        zip = "data/processed/{name}/bigscape/{name}_bigscape_as{version}.zip",
        as_dir = 'data/interim/bgcs/{name}/{version}',
        df_genomes_path = 'data/processed/{name}/tables/df_antismash_{version}_summary.csv'
    output:
        output_dir = directory("data/processed/{name}/bigscape/for_cytoscape_as-{version}")
    conda:
        "../envs/bgc_analytics.yaml"
    log: "workflow/report/logs/bigscape/bigscape_to_cytoscape/bigscape_to_cytoscape-{name}-{version}.log"
    params:
        bigscape_directory = "data/interim/bigscape/{name}_antismash_{version}/"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_bigscape_to_cytoscape.py {params.bigscape_directory} {input.as_dir} {input.df_genomes_path} {output.output_dir} 2>> {log}
        """