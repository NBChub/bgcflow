rule install_bigscape:
    output:
        bigscape = directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip)
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/ && rm BiG-SCAPE-master.zip)
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz)
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm)
        """

rule bigscape:
    input: 
        bigscape = "resources/BiG-SCAPE",
        gbk = expand("data/interim/antismash/{version}/{strains}", strains=STRAINS, version=dependency_version["antismash"])
    output:
        index = expand("data/interim/bigscape/antismash_{version}/index.html", version=dependency_version["antismash"])
    conda:
        "../envs/bigscape.yaml"
    params:
        antismash_dir = expand("data/interim/antismash/{version}", version=dependency_version["antismash"]),
        bigscape_dir = expand("data/interim/bigscape/antismash_{version}", version=dependency_version["antismash"]),
        label = "all",
    log: "workflow/report/bigscape.log"
    threads: 16
    shell:
        """
        python {input.bigscape}/bigscape.py -i {params.antismash_dir} -o {params.bigscape_dir} -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}        
        """