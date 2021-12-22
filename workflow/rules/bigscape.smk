rule install_bigscape:
    output:
        bigscape = directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/bigscape-install_bigscape.log"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip >> {log}) 
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/ && rm BiG-SCAPE-master.zip)
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz >> {log})
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm >> {log}) 
        """

rule get_bigscape_inputs:
    input:
        gbk = lambda wildcards: get_bigscape_inputs(wildcards.name, wildcards.version)
    output:
        tempdir = temp(directory("data/interim/bigscape/temp_{name}_antismash_{version}/"))
    conda:
        "../envs/bigscape.yaml"
    log: "workflow/report/logs/{name}_antismash_{version}/temp_bigscape.log"
    shell:
        """
        mkdir {output.tempdir}
        for i in {input.gbk} 
        do 
            infile=$i 
            strain=$(basename $i)
            ln -s $i "{output.tempdir}/$strain"
        done
        """

rule bigscape:
    input: 
        bigscape = "resources/BiG-SCAPE",
        antismash_dir = "data/interim/bigscape/temp_{name}_antismash_{version}/"
    output:
        index = "data/interim/bigscape/{name}_antismash_{version}/index.html"
    conda:
        "../envs/bigscape.yaml"
    params:
        bigscape_dir = "data/interim/bigscape/{name}_antismash_{version}/",
        label = "{name}_antismash_{version}",
    log: "workflow/report/logs/{name}_antismash_{version}/bigscape.log"
    threads: 16
    shell:
        """
        python {input.bigscape}/bigscape.py -i {input.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}        
        """