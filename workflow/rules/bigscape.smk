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

rule get_bigscape_inputs:
    input:
        gbk = lambda wildcards: get_antismash_inputs(wildcards.name, wildcards.version)
    output:
        tempdir = directory("data/interim/bigscape/tmp/{name}_antismash_{version}"),
        bgc_mapping = "data/interim/bigscape/{name}_antismash_{version}.csv"
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/bigscape/get_bigscape_inputs/get_bigscape_inputs-{name}_antismash_{version}.log"
    shell:
        """
        echo "Preparing BiG-SCAPE input for {wildcards.name}..." 2>> {log}
        mkdir -p {output.tempdir} 2>> {log}

        # Generate symlink for each regions in genomes in dataset
        for i in $(dirname {input.gbk})
        do
            mkdir {output.tempdir}/$(basename $i) 2>> {log}
            for r in $(ls $i/*.region*.gbk)
            do
                parent_dir=$(dirname $PWD/$r)
                filename=$(basename $r)
                (cd {output.tempdir}/$(basename $parent_dir) && ln -sf $parent_dir/$filename $filename) 2>> {log}
            done
        done

        # generate mapping for visualization
        python workflow/bgcflow/bgcflow/data/get_bigscape_mapping.py {output.tempdir} {output.bgc_mapping}
        """

rule bigscape:
    input: 
        bigscape = "resources/BiG-SCAPE",
        bgc_mapping = "data/interim/bigscape/{name}_antismash_{version}.csv"
    output:
        index = "data/interim/bigscape/{name}_antismash_{version}/index.html"
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir = "data/interim/bigscape/{name}_antismash_{version}/",
        label = "{name}_antismash_{version}",
        antismash_dir = "data/interim/bigscape/tmp/{name}_antismash_{version}/"
    log: "workflow/report/logs/bigscape/{name}_antismash_{version}/bigscape.log"
    threads: 64
    shell:
        """
        python {input.bigscape}/bigscape.py -i {params.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}
        """

rule copy_bigscape_zip:
    input: 
        bgc_mapping = "data/interim/bigscape/{name}_antismash_{version}.csv",
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