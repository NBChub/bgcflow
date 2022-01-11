rule install_bigscape:
    output:
        bigscape = directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/bigscape-install_bigscape.log"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip) > {log}
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/ && rm BiG-SCAPE-master.zip)
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz) >> {log}
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm) >> {log}
        """

rule get_bigscape_inputs:
    input:
        gbk = lambda wildcards: get_bigscape_inputs(wildcards.name, wildcards.version)
    output:
        input_list = "data/interim/bigscape/{name}_antismash_{version}.txt"
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/{name}_antismash_{version}/temp_bigscape.log"
    params:
        tempdir = directory("data/interim/bigscape/tmp/{name}_antismash_{version}")
    shell:
        """
        mkdir -p {params.tempdir} 2> {log}
        for i in {input.gbk} 
        do 
            infile=$PWD/$(dirname $i)
            strain=$(basename $infile)
            echo $infile $strain
            (cd {params.tempdir} && ln -s $infile $strain) 2>> {log}
        done
        ls {params.tempdir} > {output.input_list}
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
    log: "workflow/report/logs/{name}_antismash_{version}/bigscape.log"
    threads: 16
    shell:
        """
        python {input.bigscape}/bigscape.py -i {params.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}
        """