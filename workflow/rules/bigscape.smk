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
        input_list = "data/interim/bigscape/{name}_antismash_{version}.txt"
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
        find {output.tempdir} -iname "*.gbk" > {output.input_list}
        """

rule bigscape:
    input: 
        bigscape = "resources/BiG-SCAPE",
        input_list = "data/interim/bigscape/{name}_antismash_{version}.txt"
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